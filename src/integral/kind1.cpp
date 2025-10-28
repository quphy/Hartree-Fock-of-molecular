#include <cmath>

#include "kind1.h"

/*
kind1表示是obara_saika积分方式，详情可见JCP 84.7 (1986): 3963-3974
选择该方法进行overlap和kinetic，递推关系天然高效，适合S/T这样的少中心积分
对于核吸引和双电子部分，则递推会生成大量中间量和复杂的索引管理，同时会导致
性能下降，同时鉴于我之前的波包动力学的背景因此使用Rys方法进行这部分的处理
*/

namespace hf::inte::obara_saika {

//首先是两个单独的gto之间的overlap
double overlap_integral(const double ra[3],
                        const double rb[3],
                        const int ax, const int ay, const int az,
                        const int bx, const int by, const int bz,
                        const double alpha, const double beta) {
//如果其中一个角动量量子数小于0，那么overlap就是0（一般不会有角动量小于0，但是因为
//求导的时候可能会导致这样的结果）
if (ax < 0 || ay < 0 || az < 0 || bx < 0 || by < 0 || bz < 0) {
    return 0;
}
/*
递推关系\braket{a|b}=S_{ax,ay,az,bx,by,bz}(A,B;\alpha,\beta)，其中
\braket{a|b}为纯空间轨道积分，ax为轨道a的x方向的次方，也就是笛卡尔坐标系下
的x^{ax}，A为轨道a的中心，\alpha为轨道a的收缩系数，令p=\alpha+\beta;
P_i=\frac{\alpha·A_x+\beta·B_x}{p}，以x方向为例
S_{ax,ay,az,bx,by,bz}=(P_x-A_x)S_{ax-1,ay,az,bx,by,bz}
                     +\frac{ax-1}{2p}S_{ax-2,ay,az,bx,by,bz}
                     +\frac{bx}{2p}S_{ax-1,ay,az,bx-1,by,bz}
和
S_{ax,ay,az,bx,by,bz}=(P_x-B_x)S_{ax,ay,az,bx-1,by,bz}
                     +\frac{ax}{2p}S_{ax-1,ay,az,bx-1,by,bz}
                     +\frac{bx-1}{2p}S_{ax,ay,az,bx-2,by,bz}
y和z是一样的，最后递推一开始是
S_{0,0,0,0,0,0} = ( \frac{\pi}{\alpha+\beta} )^{3/2}
                ·exp(  -\frac{\alpha·\beta}{\alpha+\beta}
                 ·（（A_x-B_x)^2+(A_y-B_y)^2+(A_z-B_z)^2） )
*/

//开始递推
else if(ax>0){
  return (((alpha * ra[0] + beta * rb[0]) / (alpha + beta) )- ra[0]) *
           overlap_integral(ra, rb, ax-1, ay, az, bx, by, bz, alpha, beta) +
           (ax - 1) / (2.0 *(alpha + beta) )*
           overlap_integral(ra, rb, ax-2, ay, az, bx, by, bz, alpha, beta) +
           bx / (2.0 * (alpha + beta)) *
           overlap_integral(ra, rb, ax-1, ay, az, bx-1, by, bz, alpha, beta);

}

else if(ay>0){
  return (((alpha * ra[1] + beta * rb[1]) / (alpha + beta) )- ra[1]) *
           overlap_integral(ra, rb, ax, ay-1, az, bx, by, bz, alpha, beta) +
           (ay - 1) / (2.0 *(alpha + beta) )*
           overlap_integral(ra, rb, ax, ay-2, az, bx, by, bz, alpha, beta) +
           by / (2.0 * (alpha + beta)) *
           overlap_integral(ra, rb, ax, ay-1, az, bx, by-1, bz, alpha, beta);

}

else if(az>0){
  return (((alpha * ra[2] + beta * rb[2]) / (alpha + beta) )- ra[2]) *
           overlap_integral(ra, rb, ax, ay, az-1, bx, by, bz, alpha, beta) +
           (az - 1) / (2.0 *(alpha + beta) )*
           overlap_integral(ra, rb, ax, ay, az-2, bx, by, bz, alpha, beta) +
           bz / (2.0 * (alpha + beta)) *
           overlap_integral(ra, rb, ax, ay, az-1, bx, by, bz-1, alpha, beta);

}

else if(bx>0){
  return (((alpha * ra[0] + beta * rb[0]) / (alpha + beta) )- rb[0]) *
           overlap_integral(ra, rb, ax, ay, az, bx-1, by, bz, alpha, beta) +
           (bx - 1) / (2.0 *(alpha + beta) )*
           overlap_integral(ra, rb, ax, ay, az, bx-2, by, bz, alpha, beta) +
           ax / (2.0 * (alpha + beta)) *
           overlap_integral(ra, rb, ax-1, ay, az, bx-1, by, bz, alpha, beta);

}

else if(by>0){
  return (((alpha * ra[1] + beta * rb[1]) / (alpha + beta) )- rb[1]) *
           overlap_integral(ra, rb, ax, ay, az, bx, by-1, bz, alpha, beta) +
           (by - 1) / (2.0 *(alpha + beta) )*
           overlap_integral(ra, rb, ax, ay, az, bx, by-2, bz, alpha, beta) +
           ay / (2.0 * (alpha + beta)) *
           overlap_integral(ra, rb, ax, ay-1, az, bx, by-1, bz, alpha, beta);

}

else if(bz>0){
  return (((alpha * ra[2] + beta * rb[2]) / (alpha + beta) )- rb[2]) *
           overlap_integral(ra, rb, ax, ay, az, bx, by, bz-1, alpha, beta) +
           (bz - 1) / (2.0 *(alpha + beta) )*
           overlap_integral(ra, rb, ax, ay, az, bx, by, bz-2, alpha, beta) +
           az / (2.0 * (alpha + beta)) *
           overlap_integral(ra, rb, ax, ay, az-1, bx, by, bz-1, alpha, beta);

}

else{
  return sqrt(pow((M_PI/(alpha+beta)),3))*
         exp(-alpha*beta/(alpha+beta)*
         (pow(ra[0]-rb[0],2)+
         pow(ra[1]-rb[1],2)+
         pow(ra[2]-rb[2],2)) ) ;

}

}

//然后是两个inte这个命名区间内的结构体定义好的GFS的overlap
double overlap_integral(const GFS &funa, const GFS &funb){
    return overlap_integral(funa.center.memptr(),funb.center.memptr(),
                            funa.angular_number[0],funa.angular_number[1],
                            funa.angular_number[2],funb.angular_number[0],
                            funb.angular_number[1],funb.angular_number[2],
                            funa.exponent, funb.exponent)*funa.coef*funb.coef;

}

//最后获得用basis这个命名区间内的结构体包括的名字的ovelap，和h的接口一致
//也能够实现输入一个basis的名字得到overlap矩阵的目的
arma::mat overlap_integral(const basis::basis &basis){
  //创立一个overlap方阵
  arma::mat overlap(basis.n_function(),basis.n_function());

  //遍历方阵中的元素,只算上对角，因为overlap是实对称矩阵
  for (int i=0;i<basis.n_function();i++){
    for(int j=i;j<basis.n_function();j++){
      
      //basis中的每个基函数都是多个gto contracted而成，我们上面对于gto的实现都是基于
      //primitive gto进行的，因此需要再次进行展开，然后累加所有primitive gto的贡献
      const auto & function_i = basis.function[i];
      const auto n_i = function_i.coefficients.n_elem;
      const auto & function_j = basis.function[j];
      const auto n_j = function_j.coefficients.n_elem;

      double Sij=0;
      //这里的二重循环就是加和所有的primitive gto的贡献
      for (arma::uword gto_i=0;gto_i<n_i;gto_i++){
         for (arma::uword gto_j=0;gto_j<n_j;gto_j++){

           GFS funi;
           funi.center =function_i.center;
           funi.angular_number= function_i.angular_number;
           funi.exponent =function_i.exponents(gto_i);
           funi.coef =function_i.coefficients(gto_i);
           GFS funj;
           funj.center=function_j.center;
           funj.angular_number=function_j.angular_number;
           funj.exponent=function_j.exponents(gto_j);
           funj.coef=function_j.coefficients(gto_j);
          Sij+=overlap_integral(funi,funj);

        }
      }
      overlap(i,j)=Sij;
      overlap(j,i)=Sij;
    }
  }
  return overlap;
}

//动能积分和overlap其实本质一致，通过laplace求导后计算两个gto之间的overlap
arma::mat kinetic_integral(const basis::basis&basis){
  //和overlap一致
  arma::mat kinetic(basis.n_function(),basis.n_function());
  
  for(int i=0;i<basis.n_function();i++){
    for(int j=i;j<basis.n_function();j++){
      const auto & function_i = basis.function[i];
      const auto n_i = function_i.coefficients.n_elem;
      const auto & function_j = basis.function[j];
      const auto n_j = function_j.coefficients.n_elem;

      double Tij=0;
      for (arma::uword gto_i=0;gto_i<n_i;gto_i++){
         for (arma::uword gto_j=0;gto_j<n_j;gto_j++){
           GFS funi;
           funi.center =function_i.center;
           funi.angular_number= function_i.angular_number;
           funi.exponent =function_i.exponents(gto_i);
           funi.coef =function_i.coefficients(gto_i);
           GFS funj;
           funj.center=function_j.center;
           funj.angular_number=function_j.angular_number;
           funj.exponent=function_j.exponents(gto_j);
           funj.coef=function_j.coefficients(gto_j);
          //将\nabla^2作用于funj，得到新的gto
          const auto laplace_j = funj.laplace();
 
          for (const auto & gto_k: laplace_j) {
            Tij += overlap_integral(funi, gto_k);
          }
        }
      }
     kinetic(i,j)=-0.5*Tij;
     kinetic(j,i)=-0.5*Tij;
    }
  }
  return kinetic;
}
}