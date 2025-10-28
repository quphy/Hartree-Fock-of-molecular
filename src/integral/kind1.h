#ifndef obara_saika_integral
#define obara_saika_integral

#include "integral.hpp"
#include "../basis/basis.hpp"
/*
kind1表示是obara_saika积分方式，详情可见JCP 84.7 (1986): 3963-3974
选择该方法进行overlap和kinetic，递推关系天然高效，适合S/T这样的少中心积分
对于核吸引和双电子部分，则递推会生成大量中间量和复杂的索引管理，同时会导致
性能下降，因此使用Rys方法进行这部分的处理
*/
namespace hf::inte::obara_saika {
arma::mat overlap_integral(const basis::basis & basis);

arma::mat kinetic_integral(const basis::basis & basis);



}

#endif