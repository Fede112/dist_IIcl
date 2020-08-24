/*******************************************************************************
* Reads dist_Icl binary input and outputs it to terminal
* 
* The output consists of 3 columns: 
* queryID*100 + center, searchID (s), search start (ss), search end (se)
******************************************************************************/

#include <fstream>
#include <iostream>
#include <unistd.h> // getopt
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


// print a ClusteredAlignment TO STDOUT
inline void printSL( const SequenceLabel & sl) 
{
    std::cout <<  sl.sID << " " << sl.sstart << " " << sl.send <<  " " << sl.label << std::endl;   
}

int main(int argc, char** argv)
{

    //-------------------------------------------------------------------------
    // Argument parser
    //-------------------------------------------------------------------------

    int opt;
    std::string inFilename;

    while ((opt = getopt(argc, argv, "ho:")) != -1) 
    {
        switch (opt) 
        {
        // case 'o':
        //     outFilename = optarg;
        //     break;
        case 'h':
            std::cout << "Reads dist_Icl binary input file and outputs it to terminal." << std::endl;
            std::cout << "Usage: " << argv[0] << " file.bin" << std::endl;
            break;

        default: /* '?' */
            std::cout << "Usage: " << argv[0] << " file.bin" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    if (optind != argc - 1) 
    {
        std::cerr << "Expected single argument after options." << std::endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        inFilename = argv[optind];
        // std::cout << "Input: " << inFilename << std::endl;
    }

    //-------------------------------------------------------------------------

    
    std::ifstream infile (inFilename, std::ifstream::binary);
    unsigned long int length = 0;
    char * buffer = NULL;
    SequenceLabel * sca = NULL;

    // read file
    if (infile) 
    {
        // get length of file:
        infile.seekg (0, infile.end);
        length = infile.tellg();
        infile.seekg (0, infile.beg);

        buffer = new char [length];
        
        std::cerr << "Reading " << length << " characters... ";
        // read data as a block:
        infile.read (buffer,length);
        std::cerr << "all characters read successfully. \n";
        
        sca = (SequenceLabel*) buffer;
    }
    else
    {
      std::cout << "error: only " << infile.gcount() << " could be read";
    }
    infile.close();
    
    unsigned long int lines = length/sizeof(SequenceLabel);

    for (unsigned long int i = 0; i < lines; ++i)
    {
        printSL(sca[i]);
    }
    

    delete[] buffer;

  return 0;
}

