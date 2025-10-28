#include "rhf.h"
#include "rohf.h"
#include "../scf/scf.h"
#include "../scf/occ.h"
#include "../integral/kind1.h"
#include "../integral/kind2.h"
#include "../scf/mixing.h"

namespace hf::hf_cpp{
arma::mat core_hamiltonian(const geometry::atoms & atoms,
                           const basis::basis &basis){
const arma::mat kinetic =inte::obara_saika::kinetic_integral(basis);
const arma::mat nuclear_attraction=inte::rys::nuclear_attraction_integral(atoms,basis);
return kinetic+nuclear_attraction;
}

scf::FockBuilder<double>
generate_fock_builder(const basis::basis & basis,
                      const arma::mat& one_electron_integral,
                      const arma::mat& two_electron_integral ){
 const auto nao=basis.n_function();
 const arma::mat overlap=inte::obara_saika::overlap_integral(basis);
 const arma::mat &H0=one_electron_integral;

  return [nao, H0,two_electron_integral](const scf::DensityMatrix<double> & density) {
    scf::FockMatrix<double> result(arma::size(density));
    result.slice(0) =H0+arma::reshape(two_electron_integral*arma::vectorise(density.slice(0)),nao, nao);
    return result;
  };
}

scf::EnergyBuilder<double>
generate_energy_builder(const geometry::atoms & atoms,
                        const arma::mat & one_electron_integral,
                        const arma::mat & two_electron_integral) {

  const arma::mat & H0=one_electron_integral;
  const arma::vec norm_squared=arma::sum(arma::square(atoms.xyz)).t();
  arma::mat distance_squared=-2.0*atoms.xyz.t() * atoms.xyz;
  distance_squared.each_col()+=norm_squared;
  distance_squared.each_row()+=norm_squared.t();
  arma::mat inverse_r=1.0/arma::sqrt(distance_squared);
  inverse_r.diag().zeros();
  const double repulsion_among_cores=0.5*arma::dot(atoms.atomnumbers,inverse_r*atoms.atomnumbers);
  return [H0, two_electron_integral, repulsion_among_cores](
    const scf::DensityMatrix<double> & density) {
    double result=0;
    for (arma::uword i = 0; i < density.n_slices; i++) {
      result+=0.5 * arma::dot(arma::vectorise(density.slice(i)),two_electron_integral*
              arma::vectorise(density.slice(i)))+arma::dot(density.slice(i), H0);
    }
    return result+repulsion_among_cores;
  };
}

RHFSetup<double> setup(const nlohmann::json & input,
                       const geometry::atoms & atoms,
                       const basis::basis & basis) {
  scf::OverlapMatrix<double> overlap(basis.n_function(),basis.n_function(),1);
  overlap.slice(0) = inte::obara_saika::overlap_integral(basis);
  const arma::mat coulomb_integral=inte::rys::two_electron_integral(basis);
  arma::mat exchange_integral(arma::size(coulomb_integral));
  const auto nao = basis.n_function();
  for (int i = 0; i < nao; i++) {
    for (int j = 0; j < nao; j++) {
      for (int k = 0; k < nao; k++) {
        for (int l = 0; l < nao; l++) {
          exchange_integral(i+l*nao,k+j*nao)=coulomb_integral(i+j*nao,k+l*nao);
        }
      }
    }
  }
  const arma::mat two_electron_integral=coulomb_integral-0.5*exchange_integral;
  const arma::mat H0 = core_hamiltonian(atoms, basis);
  const double n_elec_per_orb=2.0;
  const scf::OccupationBuilder occupation_builder=scf::occupation::simple_occupation(n_elec_per_orb);
  const scf::FockBuilder<double> fock_builder=generate_fock_builder(basis, H0, two_electron_integral);
  const scf::EnergyBuilder<double> energy_builder=generate_energy_builder(atoms, H0, two_electron_integral);
  RHFSetup<double> setup;
  setup.atoms = atoms;
  setup.overlap = overlap;
  setup.basis = basis;
  setup.occupation_builder = occupation_builder;
  setup.fock_builder = fock_builder;
  setup.energy_builder = energy_builder;
  setup.H0 = H0;
  return setup;
}

scf::DensityMatrix<double>
initial_guess(const nlohmann::json & input,
              const arma::mat & H0,
              const scf::OverlapMatrix<double> overlap,
              const scf::OccupationBuilder occupation_builder,
              const int n_elec) {
const std::string initial_guess_method = input.at("initial_guess");
if (initial_guess_method == "H0") {
 scf::DensityMatrix<double> guess(arma::size(overlap));
 arma::cx_vec eigvals;
 arma::cx_mat eigvecs;
 arma::eig_pair(eigvals, eigvecs, H0, overlap.slice(0));
 assert(H0.is_hermitian() && overlap.slice(0).is_hermitian());
 const arma::vec real_eigvals=arma::real(eigvals);
 const arma::uvec sort_index=arma::sort_index(real_eigvals);
 const arma::vec occupation_vector=occupation_builder(real_eigvals(sort_index),n_elec);
 arma::mat sorted_orbitals = arma::real(eigvecs.cols(sort_index));
 const arma::rowvec normalization_constant=arma::real(arma::sum(sorted_orbitals%(overlap.slice(0)*sorted_orbitals)));
 sorted_orbitals.each_row()%=1.0/arma::sqrt(normalization_constant);
 const arma::mat density=sorted_orbitals * arma::diagmat(occupation_vector)*sorted_orbitals.t();
 guess.slice(0) = density;
 return guess;
} 
else {
throw Error("only H0 guess is implemented in this version");
}
}

RHFResult<double> rhf(const nlohmann::json & input,
                      const RHFSetup<double> & setup) {
  const int max_iter=input.at("max_iter");
  const double energy_tolerance=input.at("energy_tolerance");
  const scf::DensityMatrix<double>
  guess=initial_guess(input,setup.H0,setup.overlap,setup.occupation_builder,setup.atoms.numberelec());

  if(input.at("mixing")=="simple_mixing") {
   const double alpha=input.at("mixing_alpha");
   const auto update_method=scf::simple_mixing<scf::DensityMatrix<double>>{alpha, guess};

    const auto scf_result = scf::scf<scf::simple_mixing<scf::DensityMatrix<double>>,
                            double>(setup.energy_builder,setup.fock_builder,setup.occupation_builder,
                            update_method,setup.overlap,guess, {(double) setup.atoms.numberelec()},max_iter, energy_tolerance);

    return{setup, scf_result};
  } 
  else{
    throw Error("only simple mixing is implemented in this version");
  }
}

nlohmann::json rhf(const nlohmann::json & input,
                   const geometry::atoms & atoms,
                   const basis::basis & basis) {
  const auto scf_setup = setup(input, atoms, basis);
  const int multiplicity = input.value("multiplicity", 1);
  if(multiplicity==1){
  const auto result = rhf(input, scf_setup);
  nlohmann::json scf_json ;
  tr::put(scf_json, "basis_labels", basis.function_label);
  return scf_json;
  }
  else{
  return rohf(input, atoms, basis);
  }
}
}