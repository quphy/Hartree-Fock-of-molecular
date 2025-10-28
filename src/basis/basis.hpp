#ifndef basis_hf
#define basis_hf

#include <armadillo>
#include <string>

#include "../geometry/geometry.hpp"

namespace hf::basis{
struct GTO{
    //定义gto的中心
  arma::vec3 center;
    //定义轨道类型(0,0,0)为s，（1，0，0）为px
  arma::Col<int>::fixed<3> angular_number;
    //定义gto的指数系数和缩并系数
  arma::vec exponents;
  arma::vec coefficients;
} ;

struct basis {
 std::vector<GTO> function;
 std::vector<std::string> symbols;
 arma::uvec numbers;
 arma::uvec indices;
 std::vector <std::string> function_label;
 //这里定义结构体中的变量
 std::string basis_name;
//无参构造函数，方便定义basis类
 basis();
//而后含参构造函数，将对象atom和basis传进来
 basis(const geometry:: atoms & atom, const std::string &basis_set);

  int n_atoms() const;
  int n_function() const;

  arma::uvec on_atom(const arma::uword atom_index) const;
  std::vector<arma::uvec> on_atoms() const;
  basis sort_by_angular_momentum() const;

} ; 


}

#endif 