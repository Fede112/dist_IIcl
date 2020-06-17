#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>      // std::iota
#include <vector>

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


int main(int argc, char const *argv[])
{
	std::vector<uint32_t> gamma = {0  ,1  ,0  ,9  ,10  ,2   ,3,6   ,0,0};
	std::vector<double> minDist = {0.3,0.9,0.2,1  ,1   ,0.99,1,0.32,0,0};

	auto gammaSortedIdx = sort_indexes(gamma);
	

	auto gammaZeroIt = std::upper_bound(gammaSortedIdx.begin(), gammaSortedIdx.end(), 0, 
                        [&gamma = static_cast<const std::vector<uint32_t>&>(gamma)] // capture
                            (int a, int b){ return gamma[a] >= gamma[b]; }); // lambda
    

	uint32_t gammaIdxCut = (gammaZeroIt - gammaSortedIdx.begin());

	for (int i = 0; i < gammaIdxCut; ++i)
	{
		if(minDist[gammaSortedIdx[i]]==1)
			std::cout << "gamma[gammaSortedIdx[" << i << "]]: " << gamma[gammaSortedIdx[i]] << '\n';
	}

	return 0;


}