#ifndef geometryde
#define geometryde

#include <string>
#include <armadillo>
#include <vector>

namespace hf::geometry
{
    bool isEven(int a);
    
    //用结构体给定原子信息
    struct atoms
    {
     public:
     //原子符号 H，He等
       std::vector<std::string> atomsymbols;
    //对应的原子数
       arma::uvec atomnumbers;
    //原子坐标
       arma::mat xyz;
       int charge;
       int numberatom() const;
       int numberelec() const;

    };
    
} 



#endif