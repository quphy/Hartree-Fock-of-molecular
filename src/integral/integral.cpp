#include<armadillo>

#include"../basis/basis.hpp"
#include "integral.hpp"

namespace hf::inte{
//gto求导后新的gto
std::vector<GFS> GFS::gradient(const int index_cor) const{
    /* gto=x^lx y^ly z^lz coef exp(-exponent x^2)exp(-exponent y^2)exp(-exponent z^2)
       (为方便，中心为原点)对x求导后的为=x^(lx+1) y^ly z^lz -2*coef*exponent 
       exp(-exponent x^2)exp(-exponent y^2)exp(-exponent z^2)+
       x^(lx-1) y^ly z^lz coef*lx 
       exp(-exponent x^2)exp(-exponent y^2)exp(-exponent z^2)
       (lx>0),当lx=0的时候
       =x y^ly z^lz -2*coef*exponent 
       exp(-exponent x^2)exp(-exponent y^2)exp(-exponent z^2)
    */
    arma::Col<int>::fixed<3> angular_p=arma::Col<int>::fixed<3>(arma::fill::zeros);
    angular_p(index_cor)=1;
    //第一项
    std::vector<GFS> after_nabla={
     {center,angular_number+angular_p,exponent,-coef*2.0*exponent}
    };
    //第二项
  if(angular_number(index_cor) > 0) {
    after_nabla.push_back(
        {center, angular_number - angular_p, exponent, coef * angular_number(index_cor)}
        );
  }    
  
  return after_nabla;

}

std::vector<GFS> GFS::laplace() const{
// \laplace=\nabla\cdot\nabla 因此是对每个方向作用两遍
  std::vector<GFS> after_laplace;
  for(int index_cor=0; index_cor<3; index_cor++) {
    const auto first_order = this->gradient(index_cor);
    for(const auto & i_gto : first_order) {
      const auto second_order = i_gto.gradient(index_cor);
      after_laplace.insert(after_laplace.end(), second_order.begin(), second_order.end());
    }
  }

  return after_laplace;
}

}