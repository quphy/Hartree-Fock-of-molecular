#include "cal.hpp"
#include "somes/transfer.h"
#include "basis/basis.hpp"
#include "geometry/getinfo.hpp"
#include "hf/rhf.h"
#include "hf/uhf.h"
#include "hf/rohf.h"
namespace hf{

nlohmann::json cal(const nlohmann::json & input){
std::cout<<"Geometry information is read"<<std::endl;
const geometry::atoms atoms=geometry::resolve(input.at("geometry"));
std::cout<<"Basis setvi is read"<<std::endl;
const std::string basis_string=input.at("basis");
const basis::basis basis(atoms,basis_string);
const std::string method=input.at("method");
nlohmann::json output;
output["input"]=input;

if(method == "rhf") {
  std::cout<<"RHF is used"<<std::endl;
  output["output"]=hf_cpp::rhf(input, atoms, basis);
}
if(method == "uhf") {
  std::cout<<"UHF is used"<<std::endl;
  output["output"]=hf_cpp::uhf(input, atoms, basis);
}
if(method == "rohf") {
  std::cout<<"ROHF is used"<<std::endl;  
  output["output"]=hf_cpp::rohf(input, atoms, basis);
}
return output;
}

}
