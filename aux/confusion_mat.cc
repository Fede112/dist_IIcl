#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <string>
#include <vector>


struct pfStruct
{
    std::string PFAMID;
    uint16_t sstart;
    uint16_t send;
    uint32_t label;
    
    pfStruct(std::string pfid, uint16_t ss, uint16_t se, uint32_t l): PFAMID(pfid), sstart(ss), send(se), label(l) {}
    pfStruct(): PFAMID(), sstart(0), send (0), label(-9) {}

};

struct mcStruct
{
	std::string PFAMID;
	uint32_t sID;
    uint16_t sstart;
    uint16_t send;
    int label;
    
    mcStruct(std::string pfid, uint32_t s, uint16_t ss, uint16_t se, int l): PFAMID(pfid), sID(s), sstart(ss), send(se), label(l) {}
    mcStruct(): PFAMID(), sID(0), sstart(0), send (0), label(-9) {}

};


int main(int argc, char const *argv[])
{
	std::ifstream infilePF("OUR_Uniref50GT_sorted_numlabel.txt");
	std::ifstream infileMC("sequence-labeled_filtered_uniref_01.txt");

	std::string line;

	std::vector<pfStruct> dataPF;
	dataPF.reserve(14000000);
	std::vector<mcStruct> dataMC;
	dataMC.reserve(20000000);

	while (std::getline(infilePF, line))
	{
	    std::istringstream iss(line);
	    pfStruct tmp;
	    if (!(iss >> tmp.PFAMID >> tmp.sstart >> tmp.send >> tmp.label)) { break; } // error
	    dataPF.push_back(tmp);
	}


	while (std::getline(infileMC, line))
	{
	    std::istringstream iss(line);
	    mcStruct tmp2;
	    if (!(iss >> tmp2.PFAMID >> tmp2.sID >> tmp2.sstart >> tmp2.send >> tmp2.label)) { break; } // error
	    dataMC.push_back(tmp2);
	    // process pair (a,b)
	}

	for(auto elem: dataMC)
	{
		std::cout << "elem.label: " << elem.label << '\n';
	}
	
	return 0;
}


