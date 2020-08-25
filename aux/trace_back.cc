/*******************************************************************************
* k-way merge of dist_Icl output binary files.
*
* It can run using OpenMP.
*
* Important: files need to be sorted in order for the merge to succeed!
******************************************************************************/


#include <algorithm>
#include <cstring> // memcpy
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <list>
#include <random>
#include <functional>
#include <math.h>    
#include <unistd.h> // getopt
 
// #include "datatypes.h"


#define MAX_qID 4000
// #define MAX_qID 2353198020


//-------------------------------------------------------------------------
// Data types
//-------------------------------------------------------------------------


struct SequenceLabel
{
    uint32_t sID;
    uint16_t sstart;
    uint16_t send;
    int label;
    
    SequenceLabel(uint32_t s, uint16_t ss, uint16_t se, int l): sID(s), sstart(ss), send(se), label(l) {}
    SequenceLabel(): sID(0), sstart(0), send (0), label(-9) {}

};

// needed for std::upper_bound in compute_cluster_size
struct compare_label
{
    bool operator() (const SequenceLabel & left, const SequenceLabel & right)
    {
        return left.label < right.label;
    }
    bool operator() (const SequenceLabel & left, int right)
    {
        return left.label < right;
    }
    bool operator() (int left, const SequenceLabel & right)
    {   
        return left < right.label;
    }
};


// print a ClusteredAlignment TO STDOUT
inline void printSL( const SequenceLabel & sl) 
{
    std::cout <<  sl.sID << " " << sl.sstart << " " << sl.send <<  " " << sl.label << std::endl;   
}



// Small Clustered Alignment structure (only essential data...)
struct SmallCA
{
    uint32_t qID; // qID*100 + center
    uint32_t qSize;
    uint32_t sID;
    uint16_t sstart;
    uint16_t send;
    
    SmallCA(uint32_t q, uint32_t qs, uint32_t s, uint16_t ss, uint16_t se): qID(q), qSize(qs), sID(s), sstart(ss), send(se){}
    SmallCA(): qID(0), qSize(0),sID(0), sstart(0), send (0) {}

};


// needed for std::upper_bound in compute_cluster_size
struct compare_qID
{
    bool operator() (const SmallCA & left, const SmallCA & right)
    {
        return left.qID < right.qID;
    }
    bool operator() (const SmallCA & left, uint32_t right)
    {
        return left.qID < right;
    }
    bool operator() (uint32_t left, const SmallCA & right)
    {   
        return left < right.qID;
    }
};



struct compare_sID // used by auxiliary kmerge_binary
{
    bool operator() (const SmallCA & left, const SmallCA & right)
    {
        return left.sID < right.sID;
    }
    bool operator() (const SmallCA & left, uint32_t right)
    {
        return left.sID < right;
    }
    bool operator() (uint32_t left, const SmallCA & right)
    {   
        return left < right.sID;
    }
};

// print a ClusteredAlignment TO STDOUT
inline void printSCA( const SmallCA & clusterAlign) 
{
    std::cout << clusterAlign.qID << " " << clusterAlign.qSize << " " << clusterAlign.sID << " " << clusterAlign.sstart << " " << clusterAlign.send << std::endl;   
}


//-------------------------------------------------------------------------
// Functions
//-------------------------------------------------------------------------


template <typename Iterator, typename Comparator = std::less<>>
void kmerge_pow2(std::vector<Iterator> & indexes, Comparator cmp = Comparator())
{
    
    int partitions = indexes.size() - 1;
    int merges = log2(partitions);

    // merging iterations
    for (int j = 0; j < merges; ++j)
    {
        // parallel merge
        #pragma omp parallel for
        for (int i = 0; i < partitions/2; ++i)
        {
            std::inplace_merge(indexes[2*i], indexes[2*i+1], indexes[2*i+2], cmp);
        }

        // after merge, erase the middle point indexes
        for (int i = 0; i < partitions/2; ++i)
        {
            // for small number of files , k~100, it is not a bottleneck to do std::erase.
            indexes.erase(indexes.begin() + 2*i+1 - i);
        }

        // update the number of partitions
        partitions = partitions/2;
    }

    return;
}

template <typename Iterator, typename Comparator = std::less<>>
void kmerge(std::vector<Iterator> & partIndexes, Comparator cmp = Comparator() )
{

    int partitions = partIndexes.size() - 1;
    if (partitions == 1)
        return;
    else
    {
        std::vector<Iterator> mergedPartIndexes{partIndexes[0]};
        std::vector<int> offsetBase2{0};
        int accumMask = 0;
        for (int i = 0; i < 32; i++)
        {
            int mask = 1 << i;
            if ((partitions & mask) != 0)
            {
                accumMask += mask;
                mergedPartIndexes.push_back(partIndexes[0  + accumMask ]);
                offsetBase2.push_back(accumMask);
            }
        }

        for (uint64_t i = 1; i < offsetBase2.size(); ++i)
        {
            std::vector<Iterator> localPartIndexes(partIndexes.begin() + offsetBase2[i-1], 
                                                                        partIndexes.begin() + offsetBase2[i] + 1);

            kmerge_pow2(localPartIndexes, cmp);    
        }
        
        kmerge(mergedPartIndexes, cmp);
    }
    
}


void radix_sort(unsigned char * pData, uint64_t count, uint32_t record_size, 
                            uint32_t key_size, uint32_t key_offset)
{
    typedef uint8_t k_t [key_size];
    typedef uint8_t record_t [record_size];

    // index matrix
    uint64_t mIndex[key_size][256] = {0};            
    unsigned char * pTemp = new unsigned char [count*sizeof(record_t)];
    unsigned char * pDst, * pSrc;
    uint64_t i,j,m,n;
    k_t u;
    record_t v;

    // generate histograms
    for(i = 0; i < count; i++)
    {         
        std::memcpy(&u, (pData + record_size*i + key_offset), sizeof(u));
        
        for(j = 0; j < key_size; j++)
        {
            mIndex[j][(size_t)(u[j])]++;
            // mIndex[j][(size_t)(u & 0xff)]++;
            // u >>= 8;
        }       
    }
    // convert to indices 
    for(j = 0; j < key_size; j++)
    {
        n = 0;
        for(i = 0; i < 256; i++)
        {
            m = mIndex[j][i];
            mIndex[j][i] = n;
            n += m;
        }       
    }

    // radix sort
    pDst = pTemp;                       
    pSrc = pData;
    for(j = 0; j < key_size; j++)
    {
        for(i = 0; i < count; i++)
        {
            std::memcpy(&u, (pSrc + record_size*i + key_offset), sizeof(u));
            std::memcpy(&v, (pSrc + record_size*i), sizeof(v));
            m = (size_t)(u[j]);
            // m = (size_t)(u >> (j<<3)) & 0xff;
            std::memcpy(pDst + record_size*mIndex[j][m]++, &v, sizeof(v));
        }
        std::swap(pSrc, pDst);
        
    }
    delete[] pTemp;
    return;
}



/*******************************************************************************
* Loads binary file into std::vector<T>
*
* @param filename full path to file
* @param vector std::vector to store the content of the file
******************************************************************************/
template<typename T>
int load_labels(std::string filename, std::vector<T> & vector)
{
	int value;
    std::ifstream inFile (filename);
    
    size_t bytes{0};
    size_t totalEntries{0};
    vector.reserve(MAX_qID);


	// test file open   
	if (inFile) 
	{        
	    // read the elements in the file into a vector  
	    while ( inFile >> value ) 
	    {
	        vector.emplace_back(value);
	    }
	}

   
    inFile.close();

    return 0;
}


std::vector<SequenceLabel> match_search_label(SmallCA * clusterAlign,  std::vector<int> &label, const uint64_t length)
{
    auto end = clusterAlign + length;
    auto low = clusterAlign;

    std::vector<SequenceLabel> searchLabel;
    searchLabel.reserve(length);

    while (low != end)
    {
        uint32_t val = low->qID;

        // find last occurrence
        auto high = std::upper_bound(low, end, val, compare_qID());

        // compute the difference
        auto count = high - low;

        // add cluster size to the list
        for (auto p = low; p < high; ++p)
        {
            // we use val-1 instead of val=qID because the label array is shifted by one (qID 0 is not possible)
            searchLabel.emplace_back(SequenceLabel(p->sID, p->sstart, p->send, label[val-1])); // val%7 
        }

        // move to next element in vector (not immediate next)
        low = low + count;
    }

    std::cout << "low-clusterAlign: " << low-clusterAlign << '\n';
    std::cout << "searchLabel.size(): " << searchLabel.size() << '\n';

    return searchLabel;
}


int main(int argc, char *argv[])
{

	std::string labelFilename = "../reference/labels_aftermerge_ref.txt";

	SmallCA * bufferCA;

    // label vector: cluster label for each node
    std::vector<int> label;
    std::vector<SequenceLabel> searchLabel;



    int opt;
    std::string outPath={"./"};
    std::vector<std::string> filenames;

    while ((opt = getopt(argc, argv, "ho:")) != -1) 
    {
        switch (opt) 
        {
        case 'o':
            outPath = optarg;
            break;
        case 'h':
            // go to default

        default: /* '?' */
            fprintf(stderr, "Usage: %s file1 file2 ... [-o path to output files] \n",
                    argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    if (optind >= argc) {
        fprintf(stderr, "Expected argument after options\n");
        exit(EXIT_FAILURE);
    }

    while(optind < argc)
    {
        printf("name argument = %s\n", argv[optind]);
        filenames.push_back(argv[optind]);
        ++optind;
    }
    
    std::vector<SmallCA*> indexes; 
    std::vector<std::ifstream> files; 
    std::vector<uint64_t> fileLines; 
    uint64_t bytes{0};
    uint64_t lines{0};
    uint64_t totalLines{0};
    for (const auto & filename: filenames)
    {
        std::cout << filename << std::endl;
        std::ifstream file (filename, std::ifstream::binary);
        if(file.is_open())
        {
            file.seekg (0, file.end);
            bytes = file.tellg();
            lines = bytes/sizeof(SmallCA);
            // std::cout << lines << std::endl;
            fileLines.push_back(lines);
            totalLines += lines;
            file.seekg (0, file.beg);
            files.push_back(std::move(file));
        }
        else{std::cerr << "Cannot open file: " << filename << "\n";}
    }
    
    
    bufferCA = new SmallCA [totalLines];


    indexes.push_back(bufferCA);
    auto subBuffer = bufferCA;
    for (uint64_t i = 0; i < files.size(); ++i)
    {
        files[i].read ((char*) subBuffer, sizeof(SmallCA)*fileLines[i]);
        subBuffer += fileLines[i];
        indexes.push_back(subBuffer);
    }

    kmerge(indexes, compare_sID());

    std::cout << "Check if file is sorted wrt sID: " << std::is_sorted(bufferCA, bufferCA+totalLines, compare_sID()) << std::endl;

    
    // for (int i = 0; i < totalLines; ++i)
    // {
    // 	printSCA(bufferCA[i]);
    // }

    // load_labels
    load_labels(labelFilename, label);

    std::cout << "label.size(): " << label.size() << '\n';
    for (int i = 0; i < 10; ++i)
    {
    	std::cout << label[i] << std::endl;
    }



    // sort wrt qID+center
    if (    !std::is_sorted( bufferCA, bufferCA+totalLines, compare_qID() )    )
    {   
        radix_sort((unsigned char*) bufferCA, totalLines, 16, 4, 0);
    }

    std::cout << "Check if file is sorted wrt qID: " << std::is_sorted(bufferCA, bufferCA+totalLines, compare_qID()) << std::endl;

    // match_search_label
    searchLabel = match_search_label(bufferCA, label, totalLines);

    
    delete [] bufferCA;


    // sort wrt to sID
    radix_sort((unsigned char*) searchLabel.data(), searchLabel.size(), 12, 4, 0);
    
    // sort wrt label
    radix_sort((unsigned char*) searchLabel.data(), searchLabel.size(), 12, 4, 8);

    // Sequence label sorted wrt label
    std::cout << "Check if file is sorted wrt label: " << std::is_sorted(searchLabel.begin(), searchLabel.end(), compare_label()) << std::endl;





    // split output into different files ordered by sID and metaclusters

    size_t avgLines = size_t(searchLabel.size() / filenames.size());
    std::cout << "filenames.size(): " << filenames.size() << '\n';
    std::cout << "avgLines: " << avgLines << '\n';
    std::cout << "searchLabel.size(): " << searchLabel.size() << '\n';
    auto localStartIt = searchLabel.begin();
    auto localEndIt = localStartIt;
    for (int i = 1; i <= filenames.size(); ++i)
    {
    	std::string outFilename = outPath + "sequence-labels_" + std::to_string(i) + ".bin";
    	std::cout << "outFilename: " << outFilename << '\n';
    	std::ofstream FILE(outFilename, std::ofstream::binary);
		

        if ( (searchLabel.end() - (localStartIt + avgLines)) <= 0 )
        {
            std::cout << "used end!" << '\n';
            localEndIt = searchLabel.end();
        }
        else
            localEndIt = std::upper_bound(localStartIt, searchLabel.end(), (localStartIt + avgLines)->label, compare_label());

    	auto diff = localEndIt - localStartIt;
    	std::cout << "diff: " << diff << '\n';
    	std::cout << "localEndIt->label: " << localEndIt->label << '\n';

		FILE.write((char*)&(*localStartIt), (localEndIt - localStartIt) * sizeof(SequenceLabel));
		localStartIt = localEndIt;
    	
    	FILE.close();

    }


    return 0;

}