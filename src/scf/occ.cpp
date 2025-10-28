#include "occ.h"

namespace hf::scf::occupation{

OccupationBuilder simple_occupation(const double n_per_orb){
   return [n_per_orb](const arma::vec&eigenvalues, const double n){
   assert(eigenvalues.is_sorted());
  
   arma::vec occupation_vector(arma::size(eigenvalues), arma::fill::zeros);
    for (int i = 0; i < std::floor(n / n_per_orb); i++) {
      occupation_vector(i)=n_per_orb;
    }
    occupation_vector(std::ceil(n / n_per_orb))=n-arma::sum(occupation_vector);

    return occupation_vector;
  };
}
}