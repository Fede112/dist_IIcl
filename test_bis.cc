#include <vector>
#include <iostream>
#include <stdexcept>


class my_matrix {
    std::vector<std::vector<uint32_t> >m;
    public:
    my_matrix(uint32_t x, uint32_t y) 
    {
        m.resize(x, std::vector<uint32_t>(y,0));
    }
    
    std::vector<uint32_t>& operator[](uint32_t x) {
        return m[x];
    }
};

class utriag_matrix {
    std::vector<uint32_t> buffer;
    uint32_t size;
    public:
    utriag_matrix(uint32_t dim): size(dim), buffer( (dim*(dim+1))/2 ,0 ){}
    
    uint32_t& at(uint32_t i, uint32_t j)
    {
        if(i>=j){ throw std::invalid_argument( "at(i,j) with i >= j is not an element of utriag_matrix." );}
        return buffer[size*i - (i*(i+1))/2 - i + j - 1 ];
    }
};

int main(int argc, char const *argv[])
{
    utriag_matrix m(3);

    m.at(0,1) = 2;
    m.at(0,2) = 7;
    m.at(1,2) = 9;
    m.at(1,0) = 9;
    std::cout << "m[0][1]: " << m.at(0,1) << '\n';
    std::cout << "m[0][2]: " << m.at(0,2) << '\n';
    std::cout << "m[1][2]: " << m.at(1,2) << '\n';
    return 0;
}