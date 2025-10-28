#include <string>
#include <armadillo>
#include <vector>
#include <iostream>

#include "geometry.hpp"


namespace hf::geometry {
int atoms::numberatom() const{
//检查原子数
   bool check_atom1 = atomnumbers.n_elem == xyz.n_cols;
   bool check_atom2 = atomnumbers.n_elem == atomsymbols.size();

//输出原子数
return atomnumbers.n_elem;

}
//判断奇偶
bool isEven(int a) {
    return a % 2 == 0;
}

int atoms::numberelec() const{
    return arma::sum(atomnumbers)-charge;

}

}
