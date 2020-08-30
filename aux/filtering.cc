#include <cassert>
#include <algorithm> // bounds
#include <cstring> // memcpy
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>

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

// needed for std::upper_bound in compute_cluster_size
struct compare_sID
{
    bool operator() (const SequenceLabel & left, const SequenceLabel & right)
    {
        return left.sID < right.sID;
    }
    bool operator() (const SequenceLabel & left, uint32_t right)
    {
        return left.sID < right;
    }
    bool operator() (uint32_t left, const SequenceLabel & right)
    {   
        return left < right.sID;
    }
};


// print a ClusteredAlignment TO STDOUT
inline void printSL( const SequenceLabel & sl) 
{
    std::cout <<  sl.sID << " " << sl.sstart << " " << sl.send <<  " " << sl.label << std::endl;   
}


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



//-------------------------------------------------------------------------
// ALIGNMENTS DISTANCE on the SEARCH
//-------------------------------------------------------------------------
double dist(const SequenceLabel & i, const SequenceLabel & j){
    uint16_t hi, lo;
    double inte, uni;
    uint16_t istart, iend, jstart, jend;
    istart = i.sstart; iend = i.send; jstart = j.sstart; jend = j.send;
    //calculate intersection
    inte=0.0;
    hi=std::min(iend,jend);
    lo=std::max(istart,jstart);
    if(hi>lo) inte=hi-lo;
    //calculate union
    hi=std::max(iend,jend);
    lo=std::min(istart,jstart);
    uni=hi-lo;
    return (uni-inte)/uni;
}
// < 0.2

std::vector<SequenceLabel> filtering(std::vector<SequenceLabel> & sequences)
{
    
    std::vector<SequenceLabel> sequencesFiltered;
    sequencesFiltered.reserve(sequences.size()/6);

    auto end = sequences.cend();
    auto low = sequences.cbegin();

    while (low != end)
    {   
        std::vector<SequenceLabel> sameSearchBuff;
        sameSearchBuff.reserve(100000);
        uint32_t sIDref = low->sID;

        // find last occurrence
        auto high = std::upper_bound(low, end, sIDref, compare_sID());

        if ((high-1)->label != low->label)
        {
            high = std::lower_bound(low, high, low->label, compare_label());
            high += 1;
        }
    

        assert((high-1)->label == low->label);
        assert((high-1)->sID == low->sID);

        // compute the difference
        auto count = high - low;

        if(count == 1){low = low + count; continue;}
        


        // choose sequences which appear more than one
        // and that overlaps with itself.
        // (e.g. 12314 12 41 && 12314 10 39 [sID ss se])
        for (auto itr_i = low; itr_i < high; ++itr_i)
        {
            for (auto itr_j = low; itr_j < high; ++itr_j)
            {
                if(itr_i == itr_j){continue;}

                if (dist(*itr_i, *itr_j)<0.2)
                {
                    sameSearchBuff.emplace_back(*itr_i);
                    break;
                }
            }
        }

        if (sameSearchBuff.size()==0){low = low + count; continue;}

        // calculate average sequence among overlapped seq
        SequenceLabel avgSequence;
        avgSequence.sID = sIDref;
        size_t sstartAcc{0}; // avoid avgSequence overflow
        size_t sendAcc{0};
        for (const auto & seq: sameSearchBuff)
        {
            sstartAcc += seq.sstart;
            sendAcc += seq.send;
        }

        avgSequence.sstart = sstartAcc / sameSearchBuff.size();
        avgSequence.send = sendAcc / sameSearchBuff.size();

        // std::cout << "avg sequence: ";
        // printSL(avgSequence);

        // find the sequence which best represents the average sequence (dist metric) 
        double minDist = 1;
        SequenceLabel exponentSequence;
        for (auto & seq: sameSearchBuff)
        {
            auto distToAvg = dist(avgSequence, seq);
            // std::cout << "distToAvg: " << distToAvg << '\n';
            if (minDist >= distToAvg)
            {
                exponentSequence = seq;
                minDist = distToAvg;
            }
        }

        sequencesFiltered.push_back(exponentSequence); 

        // std::cout << "exponentSequence: ";
        // printSL(exponentSequence);

        // move to next element in vector (not immediate next)
        low = low + count;

    }

    return sequencesFiltered;
}


int main(int argc, char const *argv[])
{

    std::unordered_map<uint32_t, std::string> dict;
    std::vector<SequenceLabel> sequences;
    std::vector<SequenceLabel> sequencesFiltered;
    std::string fileName = argv[1];
    load_file(fileName, sequences);

    // std::cout << "Filtering!" << std::endl;
    sequencesFiltered = filtering(sequences);


    // std::cout << "sequences.size(): " << sequences.size() << '\n';
    // std::cout << "sequencesFiltered.size(): " << sequencesFiltered.size() << '\n';


    for (size_t i = 0; i < sequencesFiltered.size(); ++i)
    {
        printSL(sequencesFiltered[i]);
    }


    return 0;
}

