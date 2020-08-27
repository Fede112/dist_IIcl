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

void filtering(std::vector<SequenceLabel> & sequences)
{
    // auto length = sequences.size();
    auto end = sequences.end();
    auto low = sequences.begin();

    while (low != end)
    {
        uint32_t val = low->sID;

        // find last occurrence
        auto high = std::upper_bound(low, end, val, compare_sID());

        // compute the difference
        auto count = high - low;

        // internal analysis
        for (auto p = low; p < high; ++p)
        {
            printSL(*p);
        }

        std::cout << "count: " << high - low << "\n\n";

        // move to next element in vector (not immediate next)
        low = low + count;
    }

    return;
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



int main(int argc, char const *argv[])
{

    uint32_t key;
    std::string value;

    std::ifstream infile("dictionary.txt");
    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        

        if (!(iss >> key >> value)) { std::cout << "FAIL TO READ ALL!\n"; break; } // error

        std::cout << "key: " << key << '\n';
        std::cout << "value: " << value << '\n';
        // process pair (a,b)
    }


    // if (file.is_open()) {
    //     std::string line;
    //     while (std::getline(file, line)) {
            
    //         printf("%s", line.c_str());
    //     }
    //     file.close();
    // }

    // // TODO: Bug when inFileDict fails!
    // while (inFileDict)
    // {

    //     if (inFileDict >> key >> value)
    //     {
    //         std::cout << '\n';
    //         std::cout << key << '\n';
    //         // unsigned long ukey = std::stoul(key,nullptr,0);
    //         std::cout << value << '\n';
    //         // fill dictionary
    //     }
    //     else
    //     {
    //         throw std::runtime_error("Invalid input file");
    //     }

    // }
    // if (inFileDict.peek() != EOF) {std::cout << "ERROR\n";}

    std::unordered_map<uint32_t, std::string> dict;
    std::vector<SequenceLabel> sequences;
    std:: string fileName = argv[1];
    load_file(fileName, sequences);

    for (int i = 0; i < 10; ++i)
    {
        printSL(sequences[i]);
    }


    std::cout << "Filtering!" << std::endl;
    filtering(sequences);


    return 0;
}

