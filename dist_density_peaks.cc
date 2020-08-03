#include <algorithm>    // std::sort, std::stable_sort
#include <assert.h>
#include <chrono>
#include <cmath>        // exp
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>      // std::iota
#include <unistd.h>     // getopt
#include <vector>


// #define MAX_cID 2353198020 // maximum value for cluster ID
#define MAX_cID 3593 // maximum value for cluster ID

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

template<typename densityC, typename MinDistC>
std::vector<uint32_t> find_peaks(densityC const& density, MinDistC const& minDistance)
{
    std::vector<uint32_t> peaksIdx;

    // sort indexes in decreasing order based on gamma
    // auto gammaSortedIdx = sort_indexes(gamma);
    // find index, within gammaSortedIdx, of the first zero element in gamma.
    // auto gammaZeroIt = std::upper_bound(gammaSortedIdx.begin(), gammaSortedIdx.end(), 0, 
    //                    [&gamma = static_cast<const std::vector<double>&>(gamma)] // capture
    //                        (int a, int b){ return gamma[a] >= gamma[b]; }); // lambda    
    // uint32_t gammaCut = (gammaZeroIt - gammaSortedIdx.begin())*.5;
    // std::cout << "gammaCut: " << gammaCut << '\n';

    for (uint32_t i = 1; i < minDistance.size(); ++i)
    {
        if (density[i] > 1 && minDistance[ i ] >= 1)
        {
            peaksIdx.emplace_back(i);
        }
    }

    // if (peaksIdx.size() == 0) throw std::runtime_error("No peaks found!");
    return peaksIdx;
}


class utriag_matrix {
    uint32_t dim;
    std::vector<double> buffer;
    public:
    utriag_matrix(uint32_t d): dim(d), buffer( (d*(d+1))/2 ,0 ){}
    
    uint32_t size() {return dim;}

    double& at(uint32_t i, uint32_t j)
    {
        if(i==j){ throw std::invalid_argument( "at(i,j) with i == j is not a valid element of utriag_matrix." );}
        if (i>j){std::swap(i,j);}
        return buffer[dim*i - (i*(i+1))/2 - i + j - 1 ];
    }
};


// template <typename T>
// void ascending(T& dFirst, T& dSecond)
// {
//     if (dFirst > dSecond)
//         std::swap(dFirst, dSecond);
// }

uint32_t find(uint32_t x, std::vector<uint32_t> & labels)
{
    // find root label
    uint32_t y = x;
    while (labels[y]!=y)
        y = labels[y];

    while (labels[x] != x)
    {
        uint32_t z = labels[x];
        labels[x] = y;
        x = z;
    }
    return y;
}


void merge(uint32_t l1, uint32_t l2, std::vector<uint32_t> & labels)
{
    labels[find(l2, labels)] = find(l1, labels);
    return;
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
    std::vector<NormalizedPair> pcDistanceMat;

    // density vector
    // std::vector<uint32_t> density(MAX_cID+1, 1); // MAX_cID+1 because cID starts from 1
    std::vector<uint32_t> density(MAX_cID+1, 1); // MAX_cID+1 because cID starts from 1
    density[0]=0; // index 0 should not be used in the dataset, it is used as a void node
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

    // population of each metacluster - needed for mc merging
    typedef std::map<uint32_t, uint32_t> CounterMap;
    CounterMap mcCounts;

    // label vector: cluster label for each node
    std::vector<int> label(MAX_cID+1, -9);

    // label vector: cluster label for each node
    std::vector<double> distToPeak(MAX_cID+1, 10000.); // initialized to a number bigger than maximum distance between points


    //-------------------------------------------------------------------------
    // Argument parser
    //-------------------------------------------------------------------------

    int opt;
    std::string outFilename={"mc_output.txt"};
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
            std::cerr << "Usage: " << argv[0] << " file1 file2... [-o output merged file] \n";
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
    
    std::cout << "\nCalculating density... \n";

    begin = std::chrono::steady_clock::now();


    // loop threw input files
    for (auto filename: filenames)
    {

        std::cout << "Loading file.. ";
        load_file(filename, pcDistanceMat);
        std::cout << filename << " entries: " << pcDistanceMat.size() << '\n';

        for (auto & entry: pcDistanceMat)
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

    std::cout << "Density elapsed time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;


    //-------------------------------------------------------------------------
    // Descending sorting of indexes wrt density 
    // (auxiliary array used later to label the points)
    //-------------------------------------------------------------------------

    densitySortedId = sort_indexes(density);
 
    //-------------------------------------------------------------------------
    // Calculate min distance to higher density points
    //-------------------------------------------------------------------------

    std::cout << "\nCalculating minimum distance... \n";

    begin = std::chrono::steady_clock::now();


    // manually assign maximum minDistance to point with highest density
    link[densitySortedId[0]] = densitySortedId[0];
    minDistance[densitySortedId[0]] = 20; 
        

    // loop threw input files
    for (auto filename: filenames)
    {

	std::cout << "Loading file.. ";        
	load_file(filename,pcDistanceMat);
	std::cout << filename << " entries: " << pcDistanceMat.size() << '\n';

        for (const auto & entry: pcDistanceMat)
        {
	
	    // don't link elements with same density
        if ( density[entry.ID1] == density[entry.ID2] ) continue;

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

    std::cout << "Min. Distance elapsed time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    
    //-------------------------------------------------------------------------
    // Calculate gamma = density*minDistance
    //-------------------------------------------------------------------------

    std::cout << "\nCalculating gamma... \n";

    begin = std::chrono::steady_clock::now();
    // overflow is discarded because 0<=distance<=1
    std::transform( density.begin(), density.end(), minDistance.begin(), gamma.begin(),
                std::multiplies<double>() ); 
    end = std::chrono::steady_clock::now();
    std::cout << "Gamma elapsed time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;


    //-------------------------------------------------------------------------
    // Find peaks 
    //-------------------------------------------------------------------------
    std::cout << "Finding peaks... ";
    begin = std::chrono::steady_clock::now();
    peaksIdx = find_peaks(gamma, minDistance);
    end = std::chrono::steady_clock::now();
    std::cout << peaksIdx.size() << " peaks found \n";
    std::cout << "Finding peaks elapsed time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    
    //-------------------------------------------------------------------------
    // Label the points
    //-------------------------------------------------------------------------
    
    std::cout << "\nAssigning labels... \n";
    
    begin = std::chrono::steady_clock::now();


    // ASSIGN BASED ON DENSITY
    // // force peaks points to be root nodes, this is, they are linked to themselves
    // for (const auto & pId: peaksIdx) { link[pId] = pId; }

    // // cluster 0 will be for nodes with no pcDistanceMat entry
    // uint32_t clusterLabel = 1;
    
    // // could avoid labeling densitySortedID == 1
    // for (const auto & sIdx: densitySortedId)
    // {
    //     // peaks were forced to be root nodes
    //     if (sIdx == link[sIdx])
    //     {
    //         label[sIdx] = clusterLabel; 
    //         ++clusterLabel;
    //         continue;
    //     }
    //     // you have the same label as your parent node 
        
    //     // TODO: add threshold to the distance (it implies another read of the pcDistanceMat)
    //     label[sIdx] = label[ link[sIdx] ];
    // }


    // ASSIGN CONSIDERING MINIMUM DISTANCE TO PEAKS

    // sort peaksIdx (not necessary)
    // std::stable_sort(peaksIdx.begin(), peaksIdx.end()); //, [](uint32_t i1, uint32_t i2) {return i1 > i2;} );

    // sort peaksIdx wrt gamma
    // std::stable_sort(peaksIdx.begin(), peaksIdx.end(),
       // [&gamma](uint32_t i1, uint32_t i2) {return gamma[i1] > gamma[i2];});

    
    // label the peaks: labels go from 0 to peaksIdx.size()-1
    for (uint32_t i = 0; i < peaksIdx.size(); ++i)
    {
        label[peaksIdx[i]] = i;
    }


    // for (auto entry: peaksIdx)
    // {
        // std::cout << "p: " << entry << '\n';
    // }
    // exit(0);


    // loop threw all datapoints and compare distance to peaks
    for (auto filename: filenames)
    {

    std::cout << "Loading file.. ";        
    load_file(filename,pcDistanceMat);
    std::cout << filename << " entries: " << pcDistanceMat.size() << '\n';

        for (auto & entry: pcDistanceMat)
        {

            // auto itr1 = std::find(peaksIdx.begin(), peaksIdx.end(), entry.ID1);
            // auto itr2 = std::find(peaksIdx.begin(), peaksIdx.end(), entry.ID2);
            auto itr1 = std::lower_bound(peaksIdx.begin(), peaksIdx.end(), entry.ID1);
            if (*itr1 != entry.ID1) itr1 = peaksIdx.end();
            auto itr2 = std::lower_bound(peaksIdx.begin(), peaksIdx.end(), entry.ID2);
            if (*itr2 != entry.ID2) itr2 = peaksIdx.end();


            // if ID1 is a peak
            if( itr1 != peaksIdx.end())
            {   

                // don't assign points far away from peak
                if (entry.distance > 0.9) continue;
            
                // if ID2 is also peak then we don't do anything
                if(itr2 != peaksIdx.end()){continue;}

                if( distToPeak[entry.ID2] > entry.distance)
                {
                    distToPeak[entry.ID2] = entry.distance;
                    label[entry.ID2] = itr1 - peaksIdx.begin();
                }
                else if(distToPeak[entry.ID2] == entry.distance)
                {
                    if (gamma[ peaksIdx[label[entry.ID2]] ] < gamma[entry.ID1])
                    {
                        distToPeak[entry.ID2] = entry.distance;
                        label[entry.ID2] =  itr1 - peaksIdx.begin();                        
                    }
                }
                else continue;
            }

            // if ID2 is a peak
            if( itr2 != peaksIdx.end() ) 
            {
                // don't assign points far away from peak
                if (entry.distance > 0.9) continue;

                if( distToPeak[entry.ID1] > entry.distance)
                {
                    distToPeak[entry.ID1] = entry.distance;
                    label[entry.ID1] = itr2 - peaksIdx.begin();
                }
                else if(distToPeak[entry.ID1] == entry.distance)
                {
                    if (gamma[ peaksIdx[label[entry.ID1]] ] < gamma[entry.ID2])
                    {
                        distToPeak[entry.ID1] = entry.distance;
                        label[entry.ID1] = itr2 - peaksIdx.begin();
                    }
                }
                else continue;  

            } 
        }

    }



    end = std::chrono::steady_clock::now();

    std::cout << "Labeling elapsed time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    //-------------------------------------------------------------------------
    // Merge metaclusters
    //-------------------------------------------------------------------------
    
    std::cout << "\nMerging metaclusters... \n";
    
    begin = std::chrono::steady_clock::now();


    // mcPopulation.resize(peaksIdx.size());    
    // std::fill(mcPopulation.begin(), mcPopulation.end(), 0);

    // Count metacluster population
    // index 0 is the null position
    for (uint32_t i = 1; i < label.size(); ++i)
    {
        mcCounts[label[i]]++;
    }

    // for (auto itr = counts.cbegin(); itr != counts.cend(); ++itr)
    // {
    //     std::cout << itr->first << '\t' << itr->second << '\n';
    // }


    // metaclusters distance matrix: we only store upper triangular part
    utriag_matrix mcDistanceMat(peaksIdx.size());
    std::vector<uint32_t> labels(peaksIdx.size());
    std::iota(labels.begin(), labels.end(), 0);

    // TODO: move up these definitions with all the vector definitions

    // loop threw all datapoints and add distance to corresponding mc pair
    for (auto filename: filenames)
    {

    std::cout << "Loading file.. ";        
    load_file(filename,pcDistanceMat);
    std::cout << filename << " entries: " << pcDistanceMat.size() << '\n';

        for (auto & entry: pcDistanceMat)
        {
            if (label[entry.ID1] < 0 || label[entry.ID2] <0) continue;
            if (label[entry.ID1] == label[entry.ID2]) continue;

            mcDistanceMat.at(label[entry.ID1], label[entry.ID2]) += entry.distance;
        }

    }

    // normalize metacluster distance
    // for (uint32_t i = 0; i < 2; ++i)
    for (uint32_t i = 0; i < mcDistanceMat.size(); ++i)
    {
        for (uint32_t j = i+1; j < mcDistanceMat.size(); ++j)
        {

            // std::cout << "mcCounts[i]: " << mcCounts[i] << '\n';
            // std::cout << "mcCounts[j]: " << mcCounts[j] << '\n';
            // std::cout << "mcDistanceMat.at(i,j): " << mcDistanceMat.at(i,j) << '\n';
            mcDistanceMat.at(i,j) = (mcDistanceMat.at(i,j)) / (mcCounts[i]*mcCounts[j]);
            // std::cout << "mcDistanceMat.at(i,j): " << mcDistanceMat.at(i,j) << '\n';

            if (mcDistanceMat.at(i,j)<0.9)
            {
                if (labels[i] == labels[j]) continue;
                // merge both clusters
                merge(i,j,labels);
             }
        }
    }


    // repaint labels
    for (size_t i = 1; i < label.size(); ++i)
    {
        if (label[i] < 0) continue;
        label[i] = find(label[i],labels);
    }

    end = std::chrono::steady_clock::now();
    std::cout << "Merging elapsed time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    
    //-------------------------------------------------------------------------
    // Output
    //-------------------------------------------------------------------------

    

    // for (int i = MAX_cID - 10; i < MAX_cID+1; ++i)
    // {
    //     std::cout << "densitySortedId[i]: " << densitySortedId[i] << '\n';
    // }




    std::ofstream outFile(outFilename);
    for (size_t i = 1; i < label.size(); ++i)
    {
        // outFile << i-1 << '\t' << density[i] << '\t' << minDistance[i] << '\t' << label[i] << '\n';
        outFile << i-1 << '\t' << density[i] << '\t' << minDistance[i] << '\n';
    }

    std::ofstream outLabels("labels_own.txt");
    for (size_t i = 1; i < label.size(); ++i)
    {
        outLabels << label[i] << '\n';
    }

    std::ofstream outDistII("distII_own.txt");
    for (size_t i = 0; i < mcDistanceMat.size(); ++i)
    {
        for (size_t j = i+1; j < mcDistanceMat.size(); ++j)
        {
            outDistII << i+1 << ' ' << j+1 << ' ' << mcDistanceMat.at(i,j) << '\n';
        }
    }    


    outFile.close();
    outLabels.close();
    outDistII.close();

    endTotal = std::chrono::steady_clock::now();

    std::cout << "\nTotal elapsed time = " << std::chrono::duration_cast<std::chrono::milliseconds>(endTotal - beginTotal).count() << "[ms]" << std::endl;

    // output delta, labels, linking

    return 0;
}
