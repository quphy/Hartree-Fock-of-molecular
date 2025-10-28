#ifndef integral
#define integral

#include<armadillo>

#include"../basis/basis.hpp"

namespace hf::inte {

struct GFS{
  arma::vec3 center;
  arma::Col<int>::fixed<3> angular_number;
  double exponent;
  double coef;
  
  GFS()=default;


  std::vector<GFS> gradient(const int index_cor) const;
  std::vector<GFS> laplace() const;

};

// 2e integral
struct ERI{
  GFS A;
  GFS B;
  GFS C;
  GFS D;
};

//nai
using GFS_Pair=std::pair<GFS,GFS> ;

}

#endif