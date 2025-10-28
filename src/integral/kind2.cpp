#include <cassert>
#include <armadillo>

extern "C" {
#include <rys_roots.h>
} 


#include "kind2.h"

/*
kind2表示是rys高斯quadrature积分方式，详见JCC 29: 2722–2736, 2008
选择该方法进行核吸引和双电子部分，将多中心积分变成一维rys roots的求和，
计算形式直接且数值稳定，而对于T和S这样的少中心积分，使用roots反而会变慢，
因此使用os方法进行T和S
*/

/*
根据拉普拉斯变换1/r=\int_0^{\infity} exp(-r^2u^2)du,可以将核吸引积分
和上电子积分转化为关于t的多项式和权重的积分，t为新定义的变量，对于u从0到无穷的情况
t从0到1，详情见JCC 29: 2722–2736, 2008推导，然后对于这个对t的积分形式进行gaussian
quadrature积分
*/

namespace hf::inte::rys{

/*
在函数空间中，x^n可以视作是一组函数空间非正交的basis，比如
f(x)=\sum_{i=0} ci x^i，那么就可以表示为
f={c0,c1,c2,......}
*/
struct ryspolynomial{
arma::vec coef;
double operator()(double t2) const;
ryspolynomial operator*(double factor) const;
ryspolynomial operator*(const ryspolynomial & Qx) const;
};
/*
算符重载，计算p(x)=\sum coef_i x^i，t^2=x的值，使用的时候
比如:
ryspolynomial p 
p.coef={1.0,2.0,3.0}
double x = 2.0
double y = p(x)   
y 的值会是 1*2^0 + 2*2^1 + 3*2^2 = 17
*/    
double ryspolynomial::operator() (const double t2) const{

 arma::vec poly(coef.n_elem);
 for (arma::uword i=0;i<coef.n_elem;i++){
   poly(i)=std::pow(t2,i);
 }

 return arma::dot(coef,poly);
}

/*
同样是算符重载，这里重载*以定义数乘法和两个多项式相乘
首先先是数字乘法，这里
ryspolynomial p;
ryspolynomial q = p * 5.0;
假设p(x)=1+2x+3x^2，p的coef={1,2,3}
那么q(x)=5+10x+15x^2,p的coef={5,10,15}
*/

ryspolynomial ryspolynomial::operator*(const double fatcor) const{
   return {coef*fatcor};
}

/*
然后是两个多项式相乘ryspolynomial p和q
p=1+2x,q=3+4x+5x^2;p={1,2},q={3,4,5}
ryspolynomial p,q,g
ryspolynomial g=p*q
g=3+10x+13x^2+10x^3; g={3,10,13,10}
*/

ryspolynomial ryspolynomial::operator*(const ryspolynomial & Qx)const{
 
   arma::vec new_coef(Qx.coef.n_elem+coef.n_elem-1,arma::fill::zeros);

   for (arma::uword i=0;i<Qx.coef.n_elem;i++){
     for (arma::uword j = 0; j < coef.n_elem; j++) {
        new_coef(i + j) += coef(j) * Qx.coef(i);
    }
   }
   
   return {new_coef};
}

/*
关于多项式的递推关系，有HRR和VRR两种，详见论文中的29和30为VRR，论文中的31为HRR
VRR等于从小的角动量递推到大的角动量，HRR为把角动量的分配重新分配，总和不变，因此
我们在计算的时候，先用HRR对于一个具体输入的a,b，先通过HRR分配成a+b,0，再通过VRR
计算到0,0和一堆多项式，双电子积分就是a,b,c,d->a+b,0,c+d,0->0,0,0,0
*/

/*
现在定义结构体integralinfo，所有信息的变量名几乎和文章中的推导一致，首先是双电子
积分的信息
*/
namespace two_e{
struct integralinfo{
  int a;
  int b;
  int c;
  int d;

  double p;
  double P;
  double q;
  double Q;

  double A;
  double B;
  double C;
  double D;

  double alpha;
  double beta;
  double gamma;
  double delta;

  ryspolynomial polynomi;

  std::vector<integralinfo> hrr_b2a() const;
  std::vector<integralinfo> hrr_d2c() const;
  std::vector<integralinfo> vrr_a() const;
  std::vector<integralinfo> vrr_c() const;
};
/*
这里首先给出I的递推关系
*/
//文章中的eq31
std::vector<integralinfo> 
integralinfo::hrr_b2a() const{
if(b<=0){
  return {};
}
else{
auto rhs_1st=*this;
auto rhs_2nd=*this;
rhs_1st.a=rhs_1st.a+1;
rhs_1st.b=rhs_1st.b-1;
rhs_2nd.b=rhs_2nd.b-1;
rhs_2nd.polynomi=rhs_2nd.polynomi*(rhs_2nd.A-rhs_2nd.B);

return {rhs_1st,rhs_2nd};
}
}

//文章中的eq32
std::vector<integralinfo> 
integralinfo::hrr_d2c() const{
if(d<=0){
  return {};
}
else{
auto rhs_1st=*this;
auto rhs_2nd=*this;
rhs_1st.c=rhs_1st.c+1;
rhs_1st.d=rhs_1st.d-1;
rhs_2nd.d=rhs_2nd.d-1;
rhs_2nd.polynomi=rhs_2nd.polynomi*(rhs_2nd.C-rhs_2nd.D);

return {rhs_1st,rhs_2nd};
}
}

//文章中eq29
std::vector<integralinfo> 
integralinfo::vrr_a() const{
if(a<=0){
    return {};
}
else if(a==1){
//文章中第二项，也就是rhs_2nd=0，注意这里的lhs不是文章中的a+1，是a
auto rhs_1st=*this;
const ryspolynomial C00={arma::vec{rhs_1st.P-rhs_1st.A,
                        -rhs_1st.q * (rhs_1st.P-rhs_1st.Q)/(rhs_1st.p+rhs_1st.q)}};
rhs_1st.a=rhs_1st.a-1;
rhs_1st.polynomi=rhs_1st.polynomi*C00; 
if(c<=0){
//第三项为0
    return{rhs_1st};
}
else{
auto rhs_3rd=*this;
const ryspolynomial B00={arma::vec{0,0.5/(rhs_3rd.p+rhs_3rd.q)}};
rhs_3rd.a=rhs_3rd.a-1;
rhs_3rd.polynomi=rhs_3rd.polynomi*B00*rhs_3rd.c;
rhs_3rd.c=rhs_3rd.c-1;
 return {rhs_1st,rhs_3rd};
}
}
else{
auto rhs_1st=*this;
const ryspolynomial C00={arma::vec{rhs_1st.P-rhs_1st.A,
                        -rhs_1st.q * (rhs_1st.P-rhs_1st.Q)/(rhs_1st.p+rhs_1st.q)}};
rhs_1st.a=rhs_1st.a-1;
rhs_1st.polynomi=rhs_1st.polynomi*C00; 

auto rhs_2nd=*this;
const ryspolynomial B10={arma::vec{0.5/rhs_2nd.p,
                         -0.5*rhs_2nd.q/rhs_2nd.p/(rhs_2nd.p+rhs_2nd.q)}};
rhs_2nd.polynomi=rhs_2nd.polynomi*B10*(rhs_2nd.a-1);
rhs_2nd.a=rhs_2nd.a-2;
if(c==0){
    return{rhs_1st,rhs_2nd};
}
else{
auto rhs_3rd=*this;
const ryspolynomial B00={arma::vec{0,0.5/(rhs_3rd.p+rhs_3rd.q)}};
rhs_3rd.a=rhs_3rd.a-1;
rhs_3rd.polynomi=rhs_3rd.polynomi*B00*rhs_3rd.c;
rhs_3rd.c=rhs_3rd.c-1;
return{rhs_1st,rhs_2nd,rhs_3rd};
}
}
}


//文章中eq30
std::vector<integralinfo> 
integralinfo::vrr_c() const{
if(c<=0){
    return {};
}
else if(c==1){
//第二项=0
auto rhs_1st=*this;
const ryspolynomial D00={arma::vec{rhs_1st.Q-rhs_1st.C,
                        -rhs_1st.p * (rhs_1st.P-rhs_1st.Q)/(rhs_1st.p+rhs_1st.q)}};
rhs_1st.c=rhs_1st.c-1;
rhs_1st.polynomi=rhs_1st.polynomi*D00; 
if(a<=0){
//第三项为0
    return{rhs_1st};
}
else{
auto rhs_2nd=*this;
const ryspolynomial B00={arma::vec{0,0.5/(rhs_2nd.p+rhs_2nd.q)}};
rhs_2nd.c=rhs_2nd.c-1;
rhs_2nd.polynomi=rhs_2nd.polynomi*B00*rhs_2nd.a;
rhs_2nd.a=rhs_2nd.a-1;
 return {rhs_1st,rhs_2nd};
}
}
else{
auto rhs_1st=*this;
const ryspolynomial D00={arma::vec{rhs_1st.Q-rhs_1st.C,
                        -rhs_1st.p * (rhs_1st.P-rhs_1st.Q)/(rhs_1st.p+rhs_1st.q)}};
rhs_1st.c=rhs_1st.c-1;
rhs_1st.polynomi=rhs_1st.polynomi*D00; 

auto rhs_3rd=*this;
const ryspolynomial B01={arma::vec{0.5/rhs_3rd.q,
                         -0.5*rhs_3rd.p/rhs_3rd.q/(rhs_3rd.p+rhs_3rd.q)}};
rhs_3rd.polynomi=rhs_3rd.polynomi*B01*(rhs_3rd.c-1);
rhs_3rd.c=rhs_3rd.c-2;
if(a<=0){
    return{rhs_1st,rhs_3rd};
}
else{
auto rhs_2nd=*this;
const ryspolynomial B00={arma::vec{0,0.5/(rhs_2nd.p+rhs_2nd.q)}};
rhs_2nd.c=rhs_2nd.c-1;
rhs_2nd.polynomi=rhs_2nd.polynomi*B00*rhs_2nd.a;
rhs_2nd.a=rhs_2nd.a-1;
return{rhs_1st,rhs_2nd,rhs_3rd};
}
}
}

/*
现在讲HRR和VRR反复递归作用于I(a,b,c,d)得到一组多项式和I(0,0,0,0)的关系
由于文章中eq20的定义关系，I(0,0,0,0)=1
 */

//现在先进行HRR反复递归
std::vector<integralinfo>
hrr(const std::vector<integralinfo> & info){

  std::vector<integralinfo> hrr_result{};

  for(const auto & i:info){
   if(i.b>0){
     const auto b_ones=i.hrr_b2a();
     const auto b_recur=hrr(b_ones);
     hrr_result.insert(hrr_result.end(),b_recur.begin(),b_recur.end());
   }
   else if(i.d>0){
     const auto d_ones=i.hrr_d2c();
     const auto d_recur=hrr(d_ones);
     hrr_result.insert(hrr_result.end(),d_recur.begin(),d_recur.end());
   }
   else{
    hrr_result.push_back(i);
   }
  }
return hrr_result;
}

//然后进行vrr递归
std::vector<integralinfo>
vrr(const std::vector<integralinfo> & info){

    std::vector<integralinfo> vrr_result{};
  
  for(const auto & i:info){
    assert(i.b==0&&i.d==0);
    if(i.a>0){
     const auto a_ones=i.vrr_a();
     const auto a_recur=vrr(a_ones);
     vrr_result.insert(vrr_result.end(),a_recur.begin(),a_recur.end());
    }
    else if(i.c>0){
     const auto c_ones=i.vrr_c();
     const auto c_recur=vrr(c_ones);
     vrr_result.insert(vrr_result.end(),c_recur.begin(),c_recur.end());
    }
    else{
     vrr_result.push_back(i);
    }  
  }
return vrr_result;
}

//然后将得到的 std::vector 变成 ryspolynomial
ryspolynomial reduce2rys(const integralinfo & info){
const auto reduced=vrr(hrr({info}));

const int total_angular=info.a+info.b+info.c+info.d;

ryspolynomial result={arma::vec(total_angular+1,arma::fill::zeros)};

for(const auto & i:reduced){
  const arma::uword n_elem=i.polynomi.coef.n_elem;
  result.coef(arma::span(0,n_elem-1))+=i.polynomi.coef;
}
return result;
}
//准备工作结束，开始计算四个gto的2e积分，由于Rys多项式的根随参数x与阶数n很敏感，
//因此为了该数值代码的稳定性，weights和roots调用libcint
double two_electron_integral(const ERI & eri_info){
   //基本量的定义
  const double p=eri_info.A.exponent+eri_info.B.exponent;
  const arma::vec3 P=(eri_info.A.exponent*eri_info.A.center+
                      eri_info.B.exponent*eri_info.B.center)/p;
  const arma::vec3 B2A=eri_info.A.center-eri_info.B.center;

  const double EAB=std::exp(-eri_info.A.exponent*eri_info.B.exponent
                   /p*arma::dot(B2A,B2A));
   
  const double q=eri_info.C.exponent+eri_info.D.exponent;
  const arma::vec3 Q=(eri_info.C.exponent*eri_info.C.center+
                      eri_info.D.exponent*eri_info.D.center)/q;
  const arma::vec3 D2C=eri_info.C.center-eri_info.D.center;

  const double ECD=std::exp(-eri_info.C.exponent*eri_info.D.exponent
                   /q*arma::dot(D2C,D2C));
  const double rho=p*q/(p+q);
  const arma::vec3 Q2P=P-Q;
  const double T =rho*arma::dot(Q2P,Q2P);

  const double factor_eri=2.0*std::pow(M_PI,2.5)/p/q/std::sqrt(p+q)*EAB*ECD;
                       

  const integralinfo Ix{eri_info.A.angular_number[0],
                         eri_info.B.angular_number[0],
                         eri_info.C.angular_number[0],
                         eri_info.D.angular_number[0],
                         p,P[0],q,Q[0],
                         eri_info.A.center[0],
                         eri_info.B.center[0],
                         eri_info.C.center[0],
                         eri_info.D.center[0],
                         eri_info.A.exponent,
                         eri_info.B.exponent,
                         eri_info.C.exponent,
                         eri_info.D.exponent,
                         {arma::vec{1.0}}};

  const integralinfo Iy{eri_info.A.angular_number[1],
                         eri_info.B.angular_number[1],
                         eri_info.C.angular_number[1],
                         eri_info.D.angular_number[1],
                         p,P[1],q,Q[1],
                         eri_info.A.center[1],
                         eri_info.B.center[1],
                         eri_info.C.center[1],
                         eri_info.D.center[1],
                         eri_info.A.exponent,
                         eri_info.B.exponent,
                         eri_info.C.exponent,
                         eri_info.D.exponent,
                         {arma::vec{1.0}}};

  const integralinfo Iz{eri_info.A.angular_number[2],
                         eri_info.B.angular_number[2],
                         eri_info.C.angular_number[2],
                         eri_info.D.angular_number[2],
                         p,P[2],q,Q[2],
                         eri_info.A.center[2],
                         eri_info.B.center[2],
                         eri_info.C.center[2],
                         eri_info.D.center[2],
                         eri_info.A.exponent,
                         eri_info.B.exponent,
                         eri_info.C.exponent,
                         eri_info.D.exponent,
                         {arma::vec{1.0}}};

  const ryspolynomial Ix_po=reduce2rys(Ix);
  const ryspolynomial Iy_po=reduce2rys(Iy);
  const ryspolynomial Iz_po=reduce2rys(Iz);

  const ryspolynomial I_po=Ix_po*Iy_po*Iz_po;

/*
rys/boys的标准一维形式是
In(T)=\int_0^1 Pn(t^2) exp(-Tt)dt,Pn(t^2)是以t^2为x的多项式
采用rys quadrature
In(T)=\int_0^1 Pn(t^2) exp(-Tt)dt=\sum_{i=1}^{N} \omega_i Pn(t_i^2)
所谓的nodes和weight 分别指的就是\omega_i和t_i，该数值积分方法的准确性可以由
数学理论验证为，当你左边的Pn的最高阶为n是，N需要取(n+1)/2，上面提到了因为考虑到
nodes和weight的准确性和稳定性，此处调用libcint库来生成
*/
const int n_rys=(I_po.coef.n_elem+1)/2;
/*
在libcint中为了构造便于求根的rys多项式，会先把区间在(0,1)的t映射成区间在
(0,\infity)的u，映射关系为
t^2=u/(1+u),u=t^2/(1-t^2)
*/
arma::vec u(n_rys);
arma::vec w(n_rys);

//调用libcint
CINTrys_roots(n_rys, T, u.memptr(), w.memptr());

//映射回t^2
const arma::vec t2=u/(u+1.0);

double sum=0;

for(int i=0;i<n_rys;i++){
   sum=sum+I_po(t2[i])*w[i];
}

return sum*eri_info.A.coef*eri_info.B.coef
          *eri_info.C.coef*eri_info.D.coef*factor_eri;

}
}
/*
ERI是一个4阶张量，但是可以把它看作一个矩阵，也就是所谓的原本的两个空间，通过直积，得到
一个新的空间，也就是所谓arbitrary basis和binary basis的关系，通过这样的方式，原本
ijkl四个维度的张量变成i*j x k*l的矩阵
*/
arma::mat two_electron_integral(const basis::basis & basis){
 const auto nao=basis.n_function();
 const auto matrix_length=nao*nao;

 arma::mat eri(matrix_length,matrix_length);
 /*
 同时利用ERI的基本置换对称性\braket{ij|kl}=\brakt{ji|kl}=\braket{ij|lk}
 =\braket{kl|ij}=\braket{ji|lk}=\braket{lk|ij}=\braket{kl|ji}
 \braket{lk|ji}，8个双电子积分对称啊，一共是i和j置换、k和l置换，ij和kl置换，
 通过限制i>=j和k>=l以及只算上三角/下三角矩阵就可以利用到这三种对称性，但是
 考虑到实际情况，我们只是实现i和j置换、k和l这样的对称性
 */
 for (int i=0;i<nao;i++){
    for(int j=0;j<=i;j++){
        for(int k=0;k<nao;k++){
            for(int l=0;l<=k;l++){
                const auto & function_i=basis.function[i];
                const auto ni=function_i.coefficients.n_elem;
                const auto & function_j=basis.function[j];
                const auto nj=function_j.coefficients.n_elem;
                const auto & function_k=basis.function[k];
                const auto nk=function_k.coefficients.n_elem;
                const auto & function_l=basis.function[l];
                const auto nl=function_l.coefficients.n_elem;

                double value=0;

                for(arma::uword gto_i=0;gto_i<ni;gto_i++){
                    for(arma::uword gto_j=0;gto_j<nj;gto_j++){
                        for(arma::uword gto_k=0;gto_k<nk;gto_k++){
                            for(arma::uword gto_l=0;gto_l<nl;gto_l++){
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

                                GFS funk;
                                funk.center=function_k.center;
                                funk.angular_number=function_k.angular_number;
                                funk.exponent=function_k.exponents(gto_k);
                                funk.coef=function_k.coefficients(gto_k);

                                GFS funl;
                                funl.center=function_l.center;
                                funl.angular_number=function_l.angular_number;
                                funl.exponent=function_l.exponents(gto_l);
                                funl.coef=function_l.coefficients(gto_l);

                                const ERI eri_info{funi,funj,funk,funl};
                                value=value+two_e::two_electron_integral(eri_info);
                            }
                        }
                    }
                }
              eri(i+nao*j,k+nao*l)=value;
              eri(j+nao*i,k+nao*l)=value;
              eri(i+nao*j,l+nao*k)=value;
              eri(j+nao*i,l+nao*k)=value; 
            }
        }
    }
 }
return eri;
}

/*
nai和eri基本完全一致
*/
namespace nai{
struct integralinfo{
int a;
int b;

double p;
double P;

double A;
double B;
double C;

double alpha;
double beta; 

ryspolynomial polynomi;

std::vector<integralinfo> hrr() const;
std::vector<integralinfo> vrr() const;
};

//eq55
std::vector<integralinfo> integralinfo::hrr() const{
if(b==0){
return{};
}
else{
auto rhs_1st=*this;
auto rhs_2nd=*this;
rhs_1st.a=rhs_1st.a+1;
rhs_1st.b=rhs_1st.b-1;
rhs_2nd.polynomi=rhs_2nd.polynomi*(rhs_2nd.A-rhs_2nd.B);
rhs_2nd.b=rhs_2nd.b-1;
return {rhs_1st,rhs_2nd};
}
}

//eq52，且b=0，因为我们会先用HRR把b干掉
std::vector<integralinfo> integralinfo::vrr()const{
assert(this->b <= 0);
if(a<=0){
    return{};
}
else if(a==1){
    auto rhs_1st=*this;
    rhs_1st.a=rhs_1st.a-1;
    const ryspolynomial R1x={arma::vec{rhs_1st.P-rhs_1st.A,-(rhs_1st.P-rhs_1st.C)}};
    rhs_1st.polynomi=rhs_1st.polynomi*R1x;
    return {rhs_1st};
}
else{
    auto rhs_1st=*this;
    rhs_1st.a=rhs_1st.a-1;
    const ryspolynomial R1x={arma::vec{rhs_1st.P-rhs_1st.A,-(rhs_1st.P-rhs_1st.C)}};
    rhs_1st.polynomi=rhs_1st.polynomi*R1x;   
    auto rhs_2nd=*this;
    const ryspolynomial R2={arma::vec{0.5/rhs_2nd.p,-0.5/rhs_2nd.p}};  
    rhs_2nd.polynomi=rhs_2nd.polynomi*R2*(rhs_2nd.a-1); 
    rhs_2nd.a=rhs_2nd.a-2;

    return {rhs_1st,rhs_2nd};
}
}

//hrr递归
std::vector<integralinfo> 
hrr(const std::vector<integralinfo> &info){
std::vector<integralinfo> hrr_result{};
for (const auto &i:info){
    if(i.b>0){
        const auto b_ones=i.hrr();
        const auto b_recur=hrr(b_ones);
        hrr_result.insert(hrr_result.end(),b_recur.begin(),b_recur.end()); 
    }
    else{
        hrr_result.push_back(i);
    }
}
return hrr_result;
}

//vrr递归
std::vector<integralinfo> 
vrr(const std::vector<integralinfo> &info){
std::vector<integralinfo> vrr_result{};
for(const auto &i :info){
    assert(i.b==0);
    if(i.a>0){
        const auto a_ones=i.vrr();
        const auto a_recur=vrr(a_ones);
        vrr_result.insert(vrr_result.end(),a_recur.begin(),a_recur.end());
    }
    else{
        vrr_result.push_back(i);
    }
}
return vrr_result;
}

//然后将得到的 std::vector 变成 ryspolynomial
ryspolynomial reduce2rys(const integralinfo & info){
 const auto reduced=vrr(hrr({info}));

 const int total_angular=info.a+info.b;
 ryspolynomial result ={arma::vec(total_angular+1,arma::fill::zeros)};

for(const auto & i:reduced){
  const arma::uword n_elem=i.polynomi.coef.n_elem;
  result.coef(arma::span(0,n_elem-1))+=i.polynomi.coef;
}
return result;
}

//两个gto的nai
double nuclear_attraction_integral(const GFS_Pair&pair,
                                   const arma::vec3& core_center,
                                   double charge){
  const double p=pair.first.exponent+pair.second.exponent;
  const arma::vec3 P=(pair.first.exponent*pair.first.center+
                      pair.second.exponent*pair.second.center)/p;
  const arma::vec3 B2A=pair.first.center-pair.second.center;

  const double EAB=std::exp(-pair.first.exponent*pair.second.exponent
                             /p*arma::dot(B2A,B2A));
  const arma::vec3 T2C=P-core_center;
  const double T=p*arma::dot(T2C,T2C);

  const double factor_nai=-charge*2.0*M_PI/p*EAB;

  const integralinfo Ix={pair.first.angular_number[0],
                         pair.second.angular_number[0],
                         p,P[0],
                         pair.first.center[0],
                         pair.second.center[0],
                         core_center[0],
                         pair.first.exponent,
                         pair.second.exponent,
                         {arma::vec{1.0}}};

  const integralinfo Iy={pair.first.angular_number[1],
                         pair.second.angular_number[1],
                         p,P[1],
                         pair.first.center[1],
                         pair.second.center[1],
                         core_center[1],
                         pair.first.exponent,
                         pair.second.exponent,
                         {arma::vec{1.0}}};

  const integralinfo Iz={pair.first.angular_number[2],
                         pair.second.angular_number[2],
                         p,P[2],
                         pair.first.center[2],
                         pair.second.center[2],
                         core_center[2],
                         pair.first.exponent,
                         pair.second.exponent,
                         {arma::vec{1.0}}};

 const ryspolynomial Ix_po=reduce2rys(Ix);
 const ryspolynomial Iy_po=reduce2rys(Iy);
 const ryspolynomial Iz_po=reduce2rys(Iz);

 const ryspolynomial I_po=Ix_po*Iy_po*Iz_po;
 const int n_rys=(I_po.coef.n_elem+1)/2;

 arma::vec u(n_rys);
 arma::vec w(n_rys);

 CINTrys_roots(n_rys,T,u.memptr(),w.memptr());

 const arma::vec t2=u/(u+1.0);
 double sum=0;
 for (int i=0;i<n_rys;i++){
    sum+=I_po(t2[i])*w[i];
 }
return sum*pair.first.coef*pair.second.coef*factor_nai;
}
}




arma::mat nuclear_attraction_integral(const geometry::atoms & atoms,
                                      const basis::basis & basis){
 arma::mat nai(basis.n_function(),basis.n_function());

 for(int i=0;i<basis.n_function();i++){
  for(int j=i;j<basis.n_function();j++){
    const auto & function_i=basis.function[i];
    const auto ni=function_i.coefficients.n_elem;
    const auto & function_j=basis.function[j];
    const auto nj=function_j.coefficients.n_elem;
    
    double value=0;
    for(arma::uword gto_i=0;gto_i<ni;gto_i++){
    for(arma::uword gto_j=0;gto_j<nj;gto_j++){
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
      for (int atomk=0;atomk<atoms.numberatom();atomk++){
        const arma::vec3 center =atoms.xyz.col(atomk);
        const double charge=atoms.atomnumbers(atomk);
        value=value+nai::nuclear_attraction_integral({funi,funj},
                                                     center,charge);
      }
  }
 }
  nai(i,j)=value;
  nai(j,i)=value;
}
}

return nai;
}
}