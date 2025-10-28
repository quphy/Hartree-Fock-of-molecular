#ifndef scf_h
#define scf_h
#include <armadillo>
#include <functional>
#include <json.hpp>

#include "../somes/transfer.h"
#include "../somes/time.h"
#include "../somes/printer.h"

namespace hf::scf{
//用cube第三维是channel
template<typename T>
using FockMatrix=arma::Cube<T>;

template<typename T>
using OverlapMatrix=arma::Cube<T>;

template<typename T>
using DensityMatrix=arma::Cube<T>;

using EigenvalueVector=arma::vec;

using OccupationVector=arma::vec;

template<typename T>
using FockBuilder=std::function<FockMatrix<T>(const DensityMatrix<T> &)>;

template<typename T>
using EnergyBuilder=std::function<double(const DensityMatrix<T> &)>;

using OccupationBuilder=std::function<OccupationVector(const arma::vec&, const double)>;

struct SimpleSCFWrapper {
  double energy;
  double diff;
};


template <typename st>
Printer<st> scf_printer=[](const st&state,const double time, 
                           const int iter) -> int {
int width=18;
int precision=8;

int total_length = 0;

total_length = 6 + width * 3;

fmt::print("SCF routine\n");
print_separator(total_length);
fmt::print("{:>{}}", " iter ", 6);
fmt::print("{:>{}}", "timer / s", width);
fmt::print("{:>{}}", "energy / Hartree", width);
fmt::print("{:>{}}", "energy change", width);
fmt::print("\n");
print_separator(total_length);

fmt::print("{:^{}}",iter,6);
print(time, width, precision,">");
print(state.energy, width, precision,">");
print(state.diff,width,precision,">");
fmt::print("\n");
    
return total_length;
};

//rhf
template<typename T>
struct scfresult{
arma::mat eigenvalues;
arma::Cube<T> orbitals;
arma::mat occupations;
DensityMatrix<T> density;
OverlapMatrix<T> overlap;
FockMatrix<T> fock;
double energy;
};

template<class st>
using Update=std::function<st(const st & state)>;

template<class Updatemethod,typename T>
scfresult<T> scf(const EnergyBuilder<T>&energy_builder,
                 const FockBuilder<T>&fock_builder,
                 const OccupationBuilder&occupation_builder,
                 const Updatemethod&update_method,
                 const OverlapMatrix<T>&overlap,
                 const DensityMatrix<T>&initial,
                 const arma::vec &n_electron,
                 const int max_iter,
                 const double tol){

 const double init_energy=energy_builder(initial);

 SimpleSCFWrapper wrapper={init_energy,0};
 const int length=scf_printer<SimpleSCFWrapper>(wrapper,0.0,0);
 //开始计时
 timecounting scf_time;
 //初始化
 Updatemethod updater=update_method;
 DensityMatrix<T> previous_state=initial;
 double previous_energy=init_energy;

 const auto nao=overlap.slice(0).n_rows;

 for(int iter=1;iter <= max_iter;iter++){
  //使用上一轮的density 构造Fock
 const FockMatrix<T> fock_mat=fock_builder(previous_state);
 //构造raw 的new state
 DensityMatrix<T> raw_new_state(nao,nao,fock_mat.n_slices);
 arma::mat eigenvalues(nao,fock_mat.n_slices);
 arma::Cube<T> orbitals(nao,nao,fock_mat.n_slices);
 arma::mat occupation_vectors(nao,fock_mat.n_slices);
 for (arma::uword i = 0; i < fock_mat.n_slices; i++) {
   arma::cx_vec eigvalues;
   arma::cx_mat eigvectors;
 //调用armadillo中广义特征值求解器
   arma::eig_pair(eigvalues,eigvectors,fock_mat.slice(i),overlap.slice(i));
   //厄米体系，只需保留real
   const arma::vec real_eigenvalues=arma::real(eigvalues);
   //排序
   const arma::uvec sort_index=arma::sort_index(real_eigenvalues);
   const arma::vec sorted_eigenvalues =real_eigenvalues(sort_index);
   arma::cx_mat sorted_eigvectors = eigvectors.cols(sort_index);
   //归一化
   const arma::rowvec norm_constant=arma::real(arma::sum(sorted_eigvectors%(overlap.slice(0) * sorted_eigvectors)));
   sorted_eigvectors.each_row()%=1.0/arma::conv_to<arma::cx_mat>::from(arma::sqrt(norm_constant));
   eigenvalues.col(i) = sorted_eigenvalues;
   orbitals.slice(i) = arma::real(sorted_eigvectors);
   //通过occupation orbital构造density
   const OccupationVector occupation=occupation_builder(sorted_eigenvalues,n_electron(i));
   occupation_vectors.col(i)=occupation;
   const arma::Mat<T> new_density =orbitals.slice(i) * arma::diagmat(occupation) * orbitals.slice(i).t();
   raw_new_state.slice(i)=new_density;
 }
 // 让updater产生mixing 函数（捕获 old_state）
    std::pair<Update<DensityMatrix<T>>, Updatemethod> renewed = updater(previous_state);
 // 把 raw_new_state作为 new_state传入mixing
    const DensityMatrix<T> mixed_state = renewed.first(raw_new_state);
 //更新updater
    updater = renewed.second;
 //用混合后的密度计算能量并判断收敛
    double new_energy = energy_builder(mixed_state);
    double diff=new_energy-previous_energy;
    wrapper= {new_energy, diff};
    scf_printer<SimpleSCFWrapper>(wrapper, scf_time.timesec(), iter);
    if (std::abs(diff)<tol) {
      print_separator(length);
      return {eigenvalues, orbitals, occupation_vectors, mixed_state, overlap,fock_mat, new_energy};
    }
    previous_state=mixed_state;
    previous_energy=new_energy;
 }
 throw Error("SCF did not converge");
}
}
#endif