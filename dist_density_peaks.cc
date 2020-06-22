#include <algorithm>    // std::sort, std::stable_sort
#include <chrono>
#include <cmath>        // exp
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>      // std::iota
#include <unistd.h>     // getopt
#include <vector>


// #define MAX_cID 2353198020 // maximum value for cluster ID
#define MAX_cID 3594 // maximum value for cluster ID

/*******************************************************************************
* Data type of the binary input file
******************************************************************************/
struct NormalizedPair
{
    uint32_t ID1;
    uint32_t ID2;
    double distance;

    NormalizedPair()=default;
    NormalizedPair(uint32_t id1, uint32_t id2, double d):  ID1(id1), ID2(id2), distance(d) {}
};

/*******************************************************************************
* Loads binary file into std::vector<T>
*
* @param filename full path to file
* @param vector std::vector to store the content of the file
******************************************************************************/
template<typename T>
int load_file(std::string filename, std::vector<T> & vector)
{
    std::ifstream inFile (filename, std::ifstream::binary);
    
    size_t bytes{0};
    size_t totalEntries{0};

    if (inFile)
    {
        inFile.seekg (0, inFile.end);
        bytes = inFile.tellg();
        inFile.seekg (0, inFile.beg);
        totalEntries = bytes/sizeof(T);

        vector.resize(totalEntries);

        // read data as a block into vector:
        inFile.read(reinterpret_cast<char*>(vector.data()), bytes); // or &buf[0] for C++98
    }
    else
    {
        throw std::ios_base::failure( "File " + filename + " could not be read." );  
    }
    inFile.close();

    return 0;
}


/*******************************************************************************
* Sort indexes based on comparing values in v. Indexes are sort in decreasing order.
* 
* It uses std::stable_sort instead of std::sort to avoid unnecessary 
* index re-orderings when v contains elements of equal values.
*
* @param v full path to file
******************************************************************************/

// TODO: add overflow error check
template <typename T>
std::vector<uint32_t> sort_indexes(const std::vector<T> &v) 
{

    // initialize original index locations
    std::vector<uint32_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    std::stable_sort(idx.begin(), idx.end(),
       [&v](T i1, T i2) {return v[i1] > v[i2];});

    return idx;
}


/*******************************************************************************
* Find all entries which have a top half value of gamma and minDistance == 1.
* 
* @param gamma vector
* @param minDistance vector
******************************************************************************/

template<typename GammaC, typename MinDistC>
std::vector<uint32_t> find_peaks(GammaC const& gamma, MinDistC const& minDistance)
{
    std::vector<uint32_t> peaksIdx;

    // sort indexes in decreasing order based on gamma
    auto gammaSortedIdx = sort_indexes(gamma);

    // find index, within gammaSortedIdx, of the first zero element in gamma.
    auto gammaZeroIt = std::upper_bound(gammaSortedIdx.begin(), gammaSortedIdx.end(), 0, 
                        [&gamma = static_cast<const std::vector<uint32_t>&>(gamma)] // capture
                            (int a, int b){ return gamma[a] >= gamma[b]; }); // lambda
    
    uint32_t gammaCut = (gammaZeroIt - gammaSortedIdx.begin())*.5;
    
    std::cout << "gammaCut: " << gammaCut << '\n';

    for (uint32_t i = 0; i < gammaCut; ++i)
    {
        if (minDistance[ gammaSortedIdx[i] ] == 1)
        {
            peaksIdx.emplace_back(gammaSortedIdx[i]);
        }
    }

    if (peaksIdx.size() == 0) throw std::runtime_error("No peaks found!");
    return peaksIdx;
}




int main(int argc, char **argv)
{
    // timers
    std::chrono::steady_clock::time_point beginTotal;
    std::chrono::steady_clock::time_point endTotal;
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;


    beginTotal = std::chrono::steady_clock::now();

        
    // Sparce distance Matrix
    std::vector<NormalizedPair> distanceMat;

    // density vector
    // std::vector<uint32_t> density(MAX_cID+1, 1); // MAX_cID+1 because cID starts from 1
    std::vector<uint32_t> density(MAX_cID+1, 1); // MAX_cID+1 because cID starts from 1
    std::vector<uint32_t> densitySortedId;

    // min distance vector
    std::vector<double> minDistance(MAX_cID+1, 0);

    // link index vector: for each node 'i' it contains the ID of the closest node with higher density 'link[i]' 
    std::vector<uint32_t> link(MAX_cID+1, 0); 
    
    // gamma vector (density*minDistance)
    std::vector<uint32_t> gammaSortedId; 
    std::vector<double> gamma(MAX_cID+1, 0);

    // vector containing the peaks ids
    std::vector<uint32_t> peaksIdx;

    // label vector: cluster label for each node
    std::vector<uint32_t> label(MAX_cID+1, 0);



    //-------------------------------------------------------------------------
    // Argument parser
    //-------------------------------------------------------------------------

    int opt;
    std::string outFilename={"labels.txt"};
    std::vector<std::string> filenames;

    while ((opt = getopt(argc, argv, "ho:")) != -1) 
    {
        switch (opt) 
        {
        case 'o':
            outFilename = optarg;
            break;
        case 'h':
            // go to default

        default: /* '?' */
            std::cerr << "Usage: " << argv[0] << " file1 file2 ... [-o output merged file] \n";
            exit(EXIT_FAILURE);
        }
    }

    if (optind >= argc) {
        std::cerr << "Expected argument after options\n";
        exit(EXIT_FAILURE);
    }

    while(optind < argc)
    {
        std::cout << "argv[optind]: " << argv[optind] << '\n';
        filenames.push_back(argv[optind]);
        ++optind;
    }

    std::cout << "outFilename: " << outFilename << '\n';

    //-------------------------------------------------------------------------
    // Calculate density
    //-------------------------------------------------------------------------
    
    std::cout << "\nCalculating density ... \n";

    begin = std::chrono::steady_clock::now();


    // loop threw input files
    for (auto filename: filenames)
    {

        std::cout << "Loading file.. ";
        load_file(filename, distanceMat);
        std::cout << filename << " entries: " << distanceMat.size() << '\n';

        for (auto & entry: distanceMat)
        {
           
            // density calculated with cut-off = 0.9
            if (entry.distance < 0.9)
            {
                density[entry.ID1] += 1;
                density[entry.ID2] += 1;
            }
        }
    }

    end = std::chrono::steady_clock::now();

    std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    //-------------------------------------------------------------------------
    // Calculate min distance to higher density points
    //-------------------------------------------------------------------------

    std::cout << "\nCalculating minimum distance ... \n";

    begin = std::chrono::steady_clock::now();

    // loop threw input files
    for (auto filename: filenames)
    {

	std::cout << "Loading file.. ";        
	load_file(filename,distanceMat);
	std::cout << filename << " entries: " << distanceMat.size() << '\n';

        for (const auto & entry: distanceMat)
        {
	    auto minDensityID = std::min(entry.ID1, entry.ID2, 
		[&density](uint32_t i1,uint32_t i2){return ( density[i1] < density[i2]);});
	    
	    auto maxDensityID = std::max(entry.ID1, entry.ID2, 
               [&density](uint32_t i1, uint32_t i2){return ( density[i1] < density[i2]);});
	    
            if (link[minDensityID] == 0)
            {
                link[minDensityID] = maxDensityID;
                minDistance[minDensityID] = entry.distance;
            }

            if (link[minDensityID] != 0)
            {
                if (entry.distance < minDistance[minDensityID])
                {
                    link[minDensityID] = maxDensityID;
                    minDistance[minDensityID] = entry.distance;    
                }
            }
            
        }
    }

    end = std::chrono::steady_clock::now();

    std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    
    //-------------------------------------------------------------------------
    // Calculate gamma = density*minDistance
    //-------------------------------------------------------------------------

    std::cout << "\nCalculating gamma ... \n";

    begin = std::chrono::steady_clock::now();
    // overflow is discarded because 0<=distance<=1
    std::transform( density.begin(), density.end(), minDistance.begin(), gamma.begin(),
                std::multiplies<double>() ); 
    end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;


    //-------------------------------------------------------------------------
    // Find peaks 
    //-------------------------------------------------------------------------
    std::cout << "Finding peaks ... \n";
    begin = std::chrono::steady_clock::now();
    peaksIdx = find_peaks(gammaSortedId, minDistance);
    end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    
    //-------------------------------------------------------------------------
    // Label the points
    //-------------------------------------------------------------------------
    
    std::cout << "\nAssigning labels ... \n";
    

    begin = std::chrono::steady_clock::now();
    // Id sorted by density
    densitySortedId = sort_indexes(density);

    // force peaks points to be root nodes, this is, they are linked to themselves
    for (const auto & pId: peaksIdx) { link[pId] = pId; }

    // cluster 0 will be for nodes with no distanceMat entry
    uint32_t clusterLabel = 0;
    
    // could avoid labeling densitySortedID == 1
    for (const auto & sIdx: densitySortedId)
    {
        // peaks were forced to be root nodes
        if (sIdx == link[sIdx])
        {
            label[sIdx] = clusterLabel; 
            ++clusterLabel;
            continue;
        }
        // you have the same label as your parent node 
        
        // TODO: add threshold to the distance (it implies another read of the distanceMat)
        label[sIdx] = label[ link[sIdx] ];
    }

    end = std::chrono::steady_clock::now();

    std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    
    
    //-------------------------------------------------------------------------
    // Output
    //-------------------------------------------------------------------------

    std::cout << "\nPeaks Idx (size= " << peaksIdx.size() << ") \n";
    for (const auto peak: peaksIdx){    std::cout << peak << '\t';  }
    std::cout << '\n';


    std::ofstream outFile(outFilename);
    for (size_t i = 1; i < label.size(); ++i)
    {
        outFile << label[i] << '\n';
    }


    endTotal = std::chrono::steady_clock::now();

    std::cout << "\nTotal elapsed time = " << std::chrono::duration_cast<std::chrono::milliseconds>(endTotal - beginTotal).count() << "[ms]" << std::endl;

    // output delta, labels, linking

    return 0;
}
