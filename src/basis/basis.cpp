#include <cstddef>
#include <string>
#include <iostream>
#include <fstream>
#include <armadillo>
#include <boost/math/special_functions/factorials.hpp>
#include <json.hpp>

#include "basis.hpp"
#include "../geometry/periodic_table.h"
#include "somes/error&warning.hpp"

namespace hf::basis
{
basis::basis(){
    function = {};
    symbols = {};
    numbers = arma::uvec{};
    function_label = {};
}

basis::basis(const geometry:: atoms & atom, const std::string &basis_set){

 const std::string basis_path = "../data/basis_set_exchange/";
 std::fstream file;
 file.open(basis_path + basis_set+ ".0.json") ;
 this->basis_name = basis_set ;
  if (!file.is_open())
  {
    throw Error("basis file was not read successfully");
  }
   else {
      //将找到的basis file 转化成一个结构化的json 
      const nlohmann::json basis_data= nlohmann::json::parse(file);
      //首先定义无参
      basis Basis;
   
      std::vector <arma::uword> atom_indices_vector;
      const int numberatom = atom.numberatom();
      int i_atom;
      //遍历原子
     for (i_atom =0; i_atom <numberatom; i_atom++){
          //确认此处的原子序数
         const int atomicnumber= atom.atomnumbers[i_atom];
          //确认原子位置，作为GTO的中心
         const arma::vec3 center= atom.xyz.col(i_atom);
            //将basis file中的element和原子序数赋予json对象
            const nlohmann::json atom_info_json = 
          basis_data["elements"][std::to_string(atomicnumber) ];
    
          //遍历原子的电子壳层，i_layer就是"electron_shells"每一个shell壳循环
            for (auto i_layer: atom_info_json["electron_shells"]){
             /* coefficients和exponents不一样，对于dunning基组一组exponents可搭配
              n组coefficients 比如H的pv6z 10个exponents，一共有6组coefficients
              构成6个cgf然后，6个cgf构成一个pv6z；但是对于pople基组，永远是一组
              exponents搭配一组coefficients（其实在劈裂价键上没有本质区别，比如3-21g,
              我大可以在价层说是一组3个exponents 然后第一组cgf的coefficints为（a,b，0），
              第二组为（0，0，c）但是down下来的basis里面不是这样写的 */

              //coefficients是二维数组，行代表cgf，列代表每个cgf对应的系数，
              //行数表明劈裂该轨道用的cgf的个数
             const auto i_layer_coefficients = i_layer["coefficients"];
              //exponets是一维数组
             const auto string_exponents= i_layer["exponents"];
              //i_layer_coefficients.size()给出二维数组的行数，也就是cgf的个数，
              //这里就是对每个cgf进行循环
               for(size_t i_cgf =0; i_cgf<i_layer_coefficients.size();i_cgf++){
                   /*这里需要考虑第二个pople基组和其他基组的不同，pople基组对于同一个
                    主量子数n，采用了同一套exponet，但是组合系数以及对应的前面的x，y，z
                   根据角量子数l进行，而非pople则对不同的l也用了不同的coefficients
                   同时和上面注释说的一样，pople是一组exponet搭配一组coefficint
                   因此这里角量子数的赋值能这么写*/
                   const auto angular_momentum =i_layer["angular_momentum"] ; 
                   const bool is_pople= angular_momentum.size()>1;
                 /* 因为c++有作用域包裹，因此希望const的同时希望根据条件复制的话
                 只能用三元运算符，有点等价于（忽略作用域的问题）
                 if (is_pople){
                   const int i_l=angular_momentum[i_cgf].get<int>(); 
                 }
                 else{
                    const int i_l=angular_momentum[0].get<int>(); 
                 }
                  */
                 const int i_l = is_pople ?angular_momentum[i_cgf].get<int>(): angular_momentum[0].get<int>();
                  std::string angular_symbol;
                 if (i_l==0){
                  angular_symbol='s';
                  }
                  else if (i_l==1){
                   angular_symbol='p';
                  }
                  else if (i_l==2){
                   angular_symbol='d';
                   }
                    else if (i_l==3){
                   angular_symbol='f';
                   }
                    else if (i_l==4){
                      angular_symbol='g';
                    }
                   else{
                      throw hf::Error("invalid l number (only support [0,4])");
                     }

                    const auto string_coefficients=i_layer_coefficients[i_cgf];
                    const size_t n_functions=string_coefficients.size();
                    if (string_exponents.size() != n_functions){
                        throw hf::Error("number of exponents not equal to the number of gtos");
                     }
                    //上面的都是字符串的，下面进行double
                   arma::vec coefficients (n_functions);
                   arma::vec exponents(n_functions);
                   for (size_t i=0;i<n_functions;i++){
                     const std::string coefficients_num=string_coefficients[i];
                     const std::string exponents_num=string_exponents[i];
                     coefficients(i)=std::stod(coefficients_num);
                     exponents(i)=std::stod(exponents_num);
                     }
                   /*gto从球谐（spherical—harmonic）出发定义更方便，通式为
                   gto=N r^l Y_{lm} exp[-\alpha(r-A)^2]，但是在实际展开
                   以及电子积分过程中，我们采用笛卡尔（cartesian）gto，通式为
                   gto=\sum_{n,j,k} T_{lm,njk} (x-x_A)^n (y-y_A)^j 
                   (z-z_A)^k exp[-\alpha(r-A)^2],其中n+k+j=l
                   笛卡尔gto有(l+1)(l+2)/2,球谐gto有(2l+1)，对于l的限定条件
                   (l+1)(l+2)/2>=(2l+1)，因此可以用笛卡尔gto展开球谐gto
                   这是一个张量与不可约表示的问题，球谐gto是SO(3)群的不可约表示
                   基底，是通过群表示论构造的。笛卡尔gto是对称张量，是通过多项式
                   限制构造的，是一个冗余的可约表示形式。因此笛卡尔gto相对于球谐
                   gto来说是一个超完备的basis，需要对此进行可约分解，以l=2为例，
                   球谐gto有5个对应Y_{l,m},l=2,m=2,1,0,-1,2，
                   笛卡尔gto有6个对应x^2,y^2,z^2,xy,yz,xz对应关系如下
                   Y_{2,0}=2z^2-x^2-y^2; Y_{2,\pm1}=xz,yz; 
                   Y_{2,\pm2}=x^2-y^2,xy; Y_{0,0}=x^2+y^2+z^2
                   因此x^2,y^2,z^2,xy,yz,xz张成的函数空间V可以可约分解为
                   V=H_{l=2} \oplus H_{l=0}*/
                  for (int ax=i_l;ax>=0;ax--){
                    const std::string ax_string="x^"+std::to_string(ax);
                    for(int ay=i_l-ax;ay>=0;ay--){
                        const std::string ay_string="y^"+std::to_string(ay);
                        const std::string az_string ="z^"+std::to_string(i_l-ax-ay);
                        const int az=i_l-ax-ay;
                        const arma::Col<int>::fixed<3> angular_number={ax,ay,az};
                        
                        //归一化系数
                        const double factor=std::sqrt(
                         boost::math::factorial<double>(ax)*
                         boost::math::factorial<double>(ay)*
                         boost::math::factorial<double>(az)/
                        (boost::math::factorial<double>(ax * 2)*
                         boost::math::factorial<double>(ay * 2)*
                         boost::math::factorial<double>(az * 2))
                        );
                        const arma::vec norm_con=
                        arma::pow(2.0*exponents/M_PI,0.75) % //%在armadillo是逐元素相乘
                        arma::pow(8.0*exponents,i_l/2.0) *factor;
                        
                        std::string label=std::to_string(i_atom+1)+
                        atom.atomsymbols[i_atom];

                        label =label+angular_symbol;
                        label =label+ax_string;
                        label =label+ay_string;
                        label =label+az_string;

                        const arma::uvec eff_coefficient=arma::find(coefficients);
                        
                        const GTO functions{center,angular_number,exponents(eff_coefficient),
                                           coefficients(eff_coefficient)
                                           % norm_con(eff_coefficient)
                        };

                        function.push_back(functions);
                        function_label.push_back(label);      
                        atom_indices_vector.push_back(i_atom);                 
                    }
                  }
                }




             }

        }
    
      numbers = atom.atomnumbers;
      symbols = atom.atomsymbols;
      indices = arma::uvec{atom_indices_vector};

    
   

    }

}

int basis::n_atoms() const{
  return numbers.n_elem;
}

int basis::n_function() const{
 return function.size();
}

arma::uvec basis::on_atom(const arma::uword atom_index) const{
  return arma::find(indices==atom_index);
}

std::vector<arma::uvec> basis::on_atoms()const{
 std::vector<arma::uvec> result(n_atoms());
 for(int i=0;i<n_atoms();i++){
  result[i]=on_atom(i);
 }
 return result;
}

basis basis::sort_by_angular_momentum() const{
  const auto n_function=function.size();
  arma::uvec angular_momentum(n_function);
  for(size_t i=0; i<n_function; i++) {
    angular_momentum(i) = arma::sum(function[i].angular_number);
  } 
  basis new_basis;
  std::vector<GTO> new_function_basis(n_function);
  std::vector<std::string> new_function_labels(n_function);

  for(arma::uword l=0; l<=arma::max(angular_momentum); l++) {
    const arma::uvec find = arma::find(angular_momentum == l);
    assert(find.is_sorted());
    for(arma::uword i=0; i<find.n_elem; i++) {
      new_basis.function.push_back(function[find(i)]);
      new_basis.function_label.push_back(function_label[find(i)]);
    }
  }
  new_basis.numbers = numbers;
  new_basis.symbols = symbols;
  return new_basis;

}


}
