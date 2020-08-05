#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h> // getopt
#include <vector>
#include <numeric>

struct NormalizedPair
{
    uint32_t ID1;
    uint32_t ID2;
    double distance;

    NormalizedPair()=default;
    NormalizedPair(uint32_t id1, uint32_t id2, double d):  ID1(id1), ID2(id2), distance(d) {}
};



enum class CSVState {
    UnquotedField,
    QuotedField,
    QuotedQuote
};

std::vector<std::string> readCSVRow(const std::string &row) {
    CSVState state = CSVState::UnquotedField;
    std::vector<std::string> fields {""};
    size_t i = 0; // index of the current field
    for (char c : row) {
        switch (state) {
            case CSVState::UnquotedField:
                switch (c) {
                    case '\t': // end of field
                              fields.push_back(""); i++;
                              break;
                    case '"': state = CSVState::QuotedField;
                              break;
                    default:  fields[i].push_back(c);
                              break; }
                break;
            case CSVState::QuotedField:
                switch (c) {
                    case '"': state = CSVState::QuotedQuote;
                              break;
                    default:  fields[i].push_back(c);
                              break; }
                break;
            case CSVState::QuotedQuote:
                switch (c) {
                    case '\t': // , after closing quote
                              fields.push_back(""); i++;
                              state = CSVState::UnquotedField;
                              break;
                    case '"': // "" -> "
                              fields[i].push_back('"');
                              state = CSVState::QuotedField;
                              break;
                    default:  // end of quote
                              state = CSVState::UnquotedField;
                              break; }
                break;
        }
    }
    return fields;
}

/// Read CSV file, Excel dialect. Accept "quoted fields ""with quotes"""
std::vector<std::vector<std::string>> readCSV(std::istream &in) {
    std::vector<std::vector<std::string>> table;
    std::string row;
    while (!in.eof()) {
        std::getline(in, row);
        if (in.bad() || in.fail()) {
            break;
        }
        auto fields = readCSVRow(row);
        table.push_back(fields);
    }
    return table;
}


int main(int argc, char **argv)
{
    
    //-------------------------------------------------------------------------
    // Argumnet parser
    //-------------------------------------------------------------------------

    int opt;
    std::string outFilename={"distance.bin"};
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
        std::cout << "name argument = " << argv[optind] << '\n';
        filenames.push_back(argv[optind]);
        ++optind;
    }


    std::vector<std::vector<std::string>> coordTable{};
    std::vector<NormalizedPair> distTable{};
    std::ifstream infile (filenames[0]);
    coordTable = readCSV(infile);

    std::vector<double> a{0, 0.5, 2, 3, 4};
    std::vector<double> b{2.3, 4, 2.2, 3, 1};
 




    for (size_t entry_i = 0; entry_i < coordTable.size(); ++entry_i)
    {
        
        for (size_t entry_j = entry_i + 1; entry_j < coordTable.size(); ++entry_j)
        {

            std::vector<double> doubleVector_i(coordTable[entry_i].size());
            std::vector<double> doubleVector_j(coordTable[entry_j].size());
            std::transform(coordTable[entry_i].begin(), coordTable[entry_i].end(), doubleVector_i.begin(), [](const std::string& val)
            {
                return std::stod(val);
            });
            std::transform(coordTable[entry_j].begin(), coordTable[entry_j].end(), doubleVector_j.begin(), [](const std::string& val)
            {
                return std::stod(val);
            });

            // euclidean squared distance
            double dx = (doubleVector_j[0]-doubleVector_i[0]);
            double dy = (doubleVector_j[1]-doubleVector_i[1]);;
            double distance = dx*dx + dy*dy;
            
            NormalizedPair distRow(entry_i+1, entry_j+1, distance);
            
            distTable.push_back(distRow);
        }
    }

    std::cout << "distTable.size(): " << distTable.size() << '\n';

    std::cout << "distTable[0].ID1: " << distTable[0].ID1 << '\n';
    std::cout << "distTable[0].ID2: " << distTable[0].ID2 << '\n';
    std::cout << "distTable[0].distance: " << distTable[0].distance << '\n';
    getchar();

    std::cout << "distTable[1].ID1: " << distTable[1].ID1 << '\n';
    std::cout << "distTable[1].ID2: " << distTable[1].ID2 << '\n';
    std::cout << "distTable[1].distance: " << distTable[1].distance << '\n';
    

    std::ofstream outFile;
    outFile.open(outFilename, std::ios::out | std::ios::binary);
    if (!distTable.empty())
        outFile.write(reinterpret_cast<char*>(&distTable[0]),
                        distTable.size() * sizeof(distTable[0]));

    return 0;
}
