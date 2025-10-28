#ifndef rys_integral
#define rys_integral

#include <armadillo>

#include "integral.hpp"
#include "../geometry/geometry.hpp"

/*
kind2表示是rys高斯quadrature积分方式，详见JCC 29: 2722–2736, 2008
选择该方法进行核吸引和双电子部分，将多中心积分变成一维rys roots的求和，
计算形式直接且数值稳定，而对于T和S这样的少中心积分，使用roots反而会变慢，
因此使用os方法进行T和S
*/
namespace hf::inte::rys{



arma::mat two_electron_integral(const basis::basis & basis);

arma::mat nuclear_attraction_integral(const geometry::atoms & atoms,
                                      const basis::basis & basis);


}


#endif