#include "rhf.h"
#include "uhf.h"
#include "../scf/scf.h"
#include "../scf/occ.h"
#include "../scf/mixing.h"
#include "../integral/kind1.h"
#include "../integral/kind2.h"

namespace hf::hf_cpp{

scf::FockBuilder<double>
generate_fock_builder_uhf(const basis::basis & basis,
                          const arma::mat & one_electron_integral,
                          const arma::mat & coulomb_integral,
                          const arma::mat & exchange_integral){
  const auto nao=basis.n_function();
  const arma::mat overlap=inte::obara_saika::overlap_integral(basis);
  const arma::mat & H0=one_electron_integral;
  return [nao,H0,coulomb_integral,exchange_integral]
  (const scf::DensityMatrix<double> & density) {
   scf::FockMatrix<double> result(arma::size(density));
    const arma::vec vec_Pa=arma::vectorise(density.slice(0));
    const arma::vec vec_Pb=arma::vectorise(density.slice(1));
    const arma::vec vec_Ptot=vec_Pa + vec_Pb;
     const arma::vec J_vec=coulomb_integral * vec_Ptot;
    const arma::vec Kalpha_vec=exchange_integral * vec_Pa;
    const arma::vec Kbeta_vec=exchange_integral * vec_Pb;

    result.slice(0)=H0+arma::reshape(J_vec-Kalpha_vec,nao,nao);
    result.slice(1)=H0+arma::reshape(J_vec-Kbeta_vec,nao,nao);

    return result;  
    };
}

scf::EnergyBuilder<double>
generate_energy_builder_uhf(const geometry::atoms & atoms,
                            const arma::mat & one_electron_integral,
                            const arma::mat & coulomb_integral,
                            const arma::mat & exchange_integral){
  const arma::mat & H0=one_electron_integral;
  const arma::vec norm_squared=arma::sum(arma::square(atoms.xyz)).t();
  arma::mat distance_squared=-2.0 * atoms.xyz.t() * atoms.xyz;
  distance_squared.each_col()+=norm_squared;
  distance_squared.each_row()+=norm_squared.t();
  arma::mat inverse_r=1.0/arma::sqrt(distance_squared);
  inverse_r.diag().zeros();
  const double repulsion_among_cores =0.5 * arma::dot(atoms.atomnumbers, inverse_r * atoms.atomnumbers);
  return [H0, coulomb_integral, exchange_integral, repulsion_among_cores]
      (const scf::DensityMatrix<double> & density){
    double energy = 0.0;
    const arma::vec vec_Pa = arma::vectorise(density.slice(0));
    const arma::vec vec_Pb = arma::vectorise(density.slice(1));
    const arma::vec vec_Ptot = vec_Pa + vec_Pb;
    energy+=arma::dot(vec_Ptot, arma::vectorise(H0));
    energy+=0.5*arma::dot(vec_Ptot, coulomb_integral * vec_Ptot);
    energy-=0.5*arma::dot(vec_Pa, exchange_integral * vec_Pa);
    energy-=0.5*arma::dot(vec_Pb, exchange_integral * vec_Pb);
    return energy + repulsion_among_cores;
      };
}

std::pair<int,int> split_electrons_by_multiplicity(int total_elec, int multiplicity) {
  int S2=multiplicity-1; 
  int n_alpha=(total_elec + S2)/2;
  int n_beta =total_elec-n_alpha;
  return {n_alpha, n_beta};
}

UHFSetup<double> setup_uhf(const nlohmann::json & input,
                           const geometry::atoms & atoms,
                           const basis::basis & basis) {

  UHFSetup<double> setup;
  const arma::mat ov = inte::obara_saika::overlap_integral(basis);
  setup.overlap = scf::OverlapMatrix<double>(basis.n_function(), basis.n_function(), 2);
  setup.overlap.slice(0) = ov;
  setup.overlap.slice(1) = ov;
  const arma::mat coulomb_integral = inte::rys::two_electron_integral(basis);
  const auto nao = basis.n_function();
  arma::mat exchange_integral(arma::size(coulomb_integral));
  for (int i = 0; i < nao; ++i) {
    for (int j = 0; j < nao; ++j) {
      for (int k = 0; k < nao; ++k) {
        for (int l = 0; l < nao; ++l) {
          exchange_integral(i+l*nao,k+j*nao)=coulomb_integral(i+j*nao,k+l*nao);
        }
      }
    }
  }
  const arma::mat H0=core_hamiltonian(atoms, basis);
  const double n_elec_per_orb=1.0;
  setup.occupation_builder=scf::occupation::simple_occupation(n_elec_per_orb);
  setup.fock_builder=generate_fock_builder_uhf(basis, H0, coulomb_integral, exchange_integral);
  setup.energy_builder=generate_energy_builder_uhf(atoms, H0, coulomb_integral, exchange_integral);
  setup.atoms=atoms;
  setup.basis=basis;
  setup.H0=H0;
  return setup;
}

scf::DensityMatrix<double> initial_guess_uhf(const nlohmann::json & input,
                                             const arma::mat & H0,
                                             const scf::OverlapMatrix<double> overlap,
                                             const scf::OccupationBuilder occupation_builder,
                                             const int n_alpha,
                                             const int n_beta) {

  const std::string initial_guess_method=input.at("initial_guess");
  if (initial_guess_method != "H0") {
    throw Error("only H0 guess is implemented in this version");
  }
  scf::DensityMatrix<double> guess(arma::size(overlap)); 
  arma::cx_vec eigvals;
  arma::cx_mat eigvecs;
  arma::eig_pair(eigvals, eigvecs, H0, overlap.slice(0));
  assert(H0.is_hermitian() && overlap.slice(0).is_hermitian());
  const arma::vec real_eigvals=arma::real(eigvals);
  const arma::uvec sort_index=arma::sort_index(real_eigvals);
  const arma::vec occ_alpha=occupation_builder(real_eigvals(sort_index), (double)n_alpha);
  const arma::vec occ_beta=occupation_builder(real_eigvals(sort_index), (double)n_beta);
  arma::mat sorted_orbitals=arma::real(eigvecs.cols(sort_index));
  const arma::rowvec normalization_constant=arma::real(arma::sum(sorted_orbitals%(overlap.slice(0)*sorted_orbitals)));
  sorted_orbitals.each_row() %= 1.0 / arma::sqrt(normalization_constant);
  std::string symmetry="keep";
  if (input.contains("symmetry")) {
  symmetry=input.at("symmetry");
  }
  if(n_alpha==n_beta && symmetry=="break"){
  const double eps=1;
  arma::arma_rng::set_seed_random(); 
    arma::mat rand_alpha=arma::randn<arma::mat>(sorted_orbitals.n_rows, sorted_orbitals.n_cols);
    arma::mat rand_beta =arma::randn<arma::mat>(sorted_orbitals.n_rows, sorted_orbitals.n_cols);
    arma::mat orb_alpha = sorted_orbitals + eps * rand_alpha;
    arma::mat orb_beta  = sorted_orbitals + eps * rand_beta;
    arma::rowvec norm_alpha = arma::real(arma::sum(orb_alpha % (overlap.slice(0) * orb_alpha)));
    arma::rowvec norm_beta  = arma::real(arma::sum(orb_beta  % (overlap.slice(0) * orb_beta)));
    orb_alpha.each_row() %= 1.0 / arma::sqrt(norm_alpha);
    orb_beta.each_row()%= 1.0 / arma::sqrt(norm_beta);
    const arma::mat density_alpha = orb_alpha * arma::diagmat(occ_alpha) * orb_alpha.t();
    const arma::mat density_beta  = orb_beta  * arma::diagmat(occ_beta)  * orb_beta.t();
    arma::mat Dalpha = 0.5 * (density_alpha + density_alpha.t());
    arma::mat Dbeta  = 0.5 * (density_beta  + density_beta.t());
    double current_na = arma::trace(Dalpha * overlap.slice(0));
    double current_nb = arma::trace(Dbeta  * overlap.slice(0));
    Dalpha *= ((double)n_alpha) / current_na;
    Dbeta  *= ((double)n_beta ) / current_nb;
    Dalpha = 0.5 * (Dalpha + Dalpha.t());
    Dbeta  = 0.5 * (Dbeta  + Dbeta.t());
    guess.slice(0) = Dalpha;
    guess.slice(1) = Dbeta;
  }
    else{
  const arma::mat density_alpha = sorted_orbitals * arma::diagmat(occ_alpha) * sorted_orbitals.t();
  const arma::mat density_beta  = sorted_orbitals * arma::diagmat(occ_beta)  * sorted_orbitals.t();
  guess.slice(0) = density_alpha;
  guess.slice(1) = density_beta;
    }
  return guess;
}

template<typename T>
nlohmann::json uhf(const UHFSetup<T> & setup) {
  static_assert(std::is_same<T, double>::value, "uhf only implemented for double in this version");
}

UHFResult<double> uhf_driver(const nlohmann::json&input, const UHFSetup<double>&setup) {
  const int max_iter = input.at("max_iter");
  const double energy_tolerance = input.at("energy_tolerance");
  const int total_elec = setup.atoms.numberelec();
  const int multiplicity = input.value("multiplicity", 1);
  const auto [n_alpha, n_beta] = split_electrons_by_multiplicity(total_elec, multiplicity);
  scf::DensityMatrix<double> guess = initial_guess_uhf(input, setup.H0, setup.overlap,
                                                       setup.occupation_builder,
                                                       n_alpha, n_beta);
  if (input.at("mixing") != "simple_mixing") {
    throw Error("only simple mixing is implemented in this version");
  }
  const double alpha=input.at("mixing_alpha");
  const auto update_method=scf::simple_mixing<scf::DensityMatrix<double>>{alpha, guess};

  arma::vec n_electrons = { (double) n_alpha, (double) n_beta };

  const auto scf_result=scf::scf<scf::simple_mixing<scf::DensityMatrix<double>>,double>
                        (setup.energy_builder,
                        setup.fock_builder,
                        setup.occupation_builder,
                        update_method,
                        setup.overlap,
                        guess,
                        n_electrons,
                        max_iter, 
                        energy_tolerance);
  UHFResult<double> res;
  res.setup = setup;
  res.scf_result = scf_result;
  return res;
}


template<>
nlohmann::json uhf<double>(const UHFSetup<double> & setup) {
  nlohmann::json scf_json;
  tr::put(scf_json, "basis_labels", setup.basis.function_label);
  return scf_json;
}


nlohmann::json uhf(const nlohmann::json & input,
                   const geometry::atoms & atoms,
                   const basis::basis & basis) {
  const auto scf_setup = setup_uhf(input, atoms, basis);
  const auto result = uhf_driver(input, scf_setup);
  nlohmann::json scf_json;
  tr::put(scf_json, "basis_labels", basis.function_label);
  tr::put(scf_json, "energy", result.scf_result.energy);
  return scf_json;
}






}