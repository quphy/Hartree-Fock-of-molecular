#include <armadillo>
#include <stdexcept>
#include <memory>
#include <cmath>

#include "rohf.h"
#include "rhf.h"
#include "../scf/scf.h"
#include "../scf/occ.h"
#include "../scf/mixing.h"
#include "../integral/kind1.h"
#include "../integral/kind2.h"

namespace hf::hf_cpp{

scf::EnergyBuilder<double>
generate_energy_builder_roothaan(const geometry::atoms & atoms,
                                 const arma::mat & H0,
                                 const arma::mat & coulomb_integral,
                                 const arma::mat & exchange_integral) {


 const arma::vec norm_squared = arma::sum(arma::square(atoms.xyz)).t();
  arma::mat distance_squared = -2.0 * atoms.xyz.t() * atoms.xyz;
  distance_squared.each_col() += norm_squared;
  distance_squared.each_row() += norm_squared.t();
  arma::mat inverse_r = 1.0 / arma::sqrt(distance_squared);
  inverse_r.diag().zeros(); 
  const double repulsion_among_cores=0.5*arma::dot(atoms.atomnumbers,inverse_r*atoms.atomnumbers);
  return [H0, coulomb_integral, exchange_integral, repulsion_among_cores]
    (const scf::DensityMatrix<double> & density) {
    const arma::mat P_alpha = density.slice(0);
    const arma::mat P_beta  = density.slice(1);
    const arma::mat P_total = P_alpha + P_beta;
    const arma::uword n_ao = P_total.n_rows;
    double energy = arma::trace(H0 * P_total);
    arma::vec vecPtot = arma::vectorise(P_total);
    arma::mat J = arma::reshape(coulomb_integral * vecPtot, n_ao, n_ao);
    double ecoul = arma::trace(P_total * J);
    arma::mat K_alpha = arma::reshape(exchange_integral * arma::vectorise(P_alpha), n_ao, n_ao);
    arma::mat K_beta  = arma::reshape(exchange_integral * arma::vectorise(P_beta), n_ao, n_ao);
    double exa = arma::trace(P_alpha * K_alpha);
    double exb = arma::trace(P_beta * K_beta);
    double eexch = -0.5 * (exa + exb);
    energy += 0.5 * ecoul + eexch;  
    return energy + repulsion_among_cores;
  };
}


scf::DensityMatrix<double>
initial_guess_rohf(const nlohmann::json & input,
                   const arma::mat & H0,
                   const scf::OverlapMatrix<double> & overlap,
                   const geometry::atoms & atoms) {
  const int multiplicity = input.value("multiplicity", 1);
  const int n_elec = atoms.numberelec();
  const int n_unpaired = multiplicity - 1;
  const int n_alpha = (n_elec + n_unpaired) / 2;
  const int n_beta  = n_elec - n_alpha;
  arma::cx_vec eigvals;
  arma::cx_mat eigvecs;
  arma::eig_pair(eigvals, eigvecs, H0, overlap.slice(0));
  const arma::vec real_eigvals = arma::real(eigvals);
  const arma::uvec sort_index = arma::sort_index(real_eigvals);
  arma::mat sorted_orbitals = arma::real(eigvecs.cols(sort_index));
  const arma::rowvec norm_const = arma::real(arma::sum(sorted_orbitals % (overlap.slice(0) * sorted_orbitals)));
  arma::rowvec tmp = norm_const;
  for (arma::uword c = 0; c < tmp.n_elem; ++c) {
    if (std::abs(tmp(c)) < 1e-14) tmp(c) = 1.0;
    sorted_orbitals.col(c) /= std::sqrt(tmp(c));
  }
  const int n_orb = (int) sorted_orbitals.n_cols;
  arma::vec occ_alpha = arma::zeros<arma::vec>(n_orb);
  arma::vec occ_beta  = arma::zeros<arma::vec>(n_orb);
  for (int i = 0; i < n_beta && i < n_orb; ++i) {
    occ_alpha(i) = 1.0;
    occ_beta(i) = 1.0;
  }
  for (int i = n_beta; i < n_alpha && i < n_orb; ++i) {
    occ_alpha(i) = 1.0;
  }
  scf::DensityMatrix<double> guess(arma::size(overlap));
  arma::mat P_alpha = sorted_orbitals * arma::diagmat(occ_alpha) * sorted_orbitals.t();
  arma::mat P_beta  = sorted_orbitals * arma::diagmat(occ_beta)  * sorted_orbitals.t();
  guess.slice(0) = P_alpha;
  guess.slice(1) = P_beta;
  return guess;
}

nlohmann::json rohf(const nlohmann::json & input,
                        const geometry::atoms & atoms,
                        const basis::basis & basis) {
  const auto n_ao = (arma::uword) basis.n_function();
  scf::OverlapMatrix<double> overlap(n_ao, n_ao, 2);
  arma::mat S = inte::obara_saika::overlap_integral(basis);
  overlap.slice(0) = S;
  overlap.slice(1) = S;
  const arma::mat coulomb_integral = inte::rys::two_electron_integral(basis);
  arma::mat exchange_integral(arma::size(coulomb_integral));
  for (int i = 0; i < (int)n_ao; i++) {
    for (int j = 0; j < (int)n_ao; j++) {
      for (int k = 0; k < (int)n_ao; k++) {
        for (int l = 0; l < (int)n_ao; l++) {
          exchange_integral(i + l * n_ao, k + j * n_ao) =
              coulomb_integral(i + j * n_ao, k + l * n_ao);
        }
      }
    }
  }
  const arma::mat H0 = core_hamiltonian(atoms, basis);
  const scf::EnergyBuilder<double> energy_builder =
      generate_energy_builder_roothaan(atoms, H0, coulomb_integral, exchange_integral);
  auto shared_mo_energy = std::make_shared<arma::vec>(arma::vec(n_ao, arma::fill::zeros));
  auto shared_mo_ea = std::make_shared<arma::vec>(arma::vec(n_ao, arma::fill::zeros));
  auto shared_mo_eb = std::make_shared<arma::vec>(arma::vec(n_ao, arma::fill::zeros));
  const int multiplicity = input.value("multiplicity", 1);
  const int n_elec_total = atoms.numberelec();
  const int n_unpaired = multiplicity - 1;
  const int n_alpha = (n_elec_total + n_unpaired) / 2;
  const int n_beta  = n_elec_total - n_alpha;
  scf::FockBuilder<double> fock_builder = [=](const scf::DensityMatrix<double> & density) mutable -> scf::FockMatrix<double> {
    arma::mat P_alpha = density.slice(0);
    arma::mat P_beta  = density.slice(1);
    arma::vec vecPtot = arma::vectorise(P_alpha + P_beta);
    arma::mat J = arma::reshape(coulomb_integral * vecPtot, n_ao, n_ao);
    arma::mat K_alpha = arma::reshape(exchange_integral * arma::vectorise(P_alpha), n_ao, n_ao);
    arma::mat K_beta  = arma::reshape(exchange_integral * arma::vectorise(P_beta),  n_ao, n_ao);
    arma::mat Fa = H0 + J - K_alpha;
    arma::mat Fb = H0 + J - K_beta;
    arma::mat Fc = 0.5 * (Fa + Fb);
    arma::mat pc = P_beta * S;
    arma::mat po = (P_alpha - P_beta) * S;
    arma::mat pv = arma::eye<arma::mat>(n_ao, n_ao) - P_alpha * S;
    arma::mat f_eff = arma::zeros<arma::mat>(n_ao, n_ao);
    f_eff += 0.5 * pc.t() * Fc * pc;
    f_eff += 0.5 * po.t() * Fc * po;
    f_eff += 0.5 * pv.t() * Fc * pv;
    f_eff += po.t() * Fb * pc;
    f_eff += po.t() * Fa * pv;
    f_eff += pv.t() * Fc * pc;
    f_eff = 0.5 * (f_eff + f_eff.t());
    arma::cx_vec eigvals;
    arma::cx_mat eigvecs;
    arma::eig_pair(eigvals, eigvecs, f_eff, S);
    arma::vec real_eigvals = arma::real(eigvals);
    arma::uvec sort_idx = arma::sort_index(real_eigvals);
    arma::vec sorted_eigs = real_eigvals(sort_idx);
    arma::cx_mat sorted_vecs = eigvecs.cols(sort_idx);
    arma::rowvec norm_const = arma::real(arma::sum(sorted_vecs % (S * sorted_vecs)));
    for (arma::uword col = 0; col < (arma::uword)norm_const.n_elem; ++col) {
      double val = norm_const(col);
      if (std::abs(val) < 1e-14) val = 1.0;
      double denom = std::sqrt(val);
      sorted_vecs.col(col) /= denom;
    }
    arma::mat C = arma::real(sorted_vecs);
    arma::mat Ca = C.t() * Fa * C;
    arma::mat Cb = C.t() * Fb * C;
    arma::vec mo_ea_local = arma::real(arma::diagvec(Ca));
    arma::vec mo_eb_local = arma::real(arma::diagvec(Cb));
    arma::vec mo_energy_local = sorted_eigs;
    *shared_mo_energy = mo_energy_local;
    *shared_mo_ea = mo_ea_local;
    *shared_mo_eb = mo_eb_local;
    scf::FockMatrix<double> result(n_ao, n_ao, 2);
    result.slice(0) = f_eff;
    result.slice(1) = f_eff;
    return result;
  };

scf::OccupationBuilder occupation_builder = [shared_mo_energy, shared_mo_ea, shared_mo_eb, n_alpha, n_beta]
  (const arma::vec & , const double spin_n_elec) -> arma::vec {
  const arma::vec & mo_e = *shared_mo_energy; 
  const arma::vec & moa  = *shared_mo_ea;     
  const int nmo = (int) mo_e.n_elem; 
  const int nocc  = std::max(n_alpha, n_beta); 
  const int ncore = std::min(n_alpha, n_beta); 
  const int nopen = nocc - ncore;             
  arma::vec occ_alpha(nmo, arma::fill::zeros);
  arma::vec occ_beta(nmo,  arma::fill::zeros);
  arma::uvec order_by_mo = arma::sort_index(mo_e);
  if (ncore > 0) {
    arma::uvec core_idx = order_by_mo.subvec(0, (arma::uword)std::max(0, ncore-1));
    for (arma::uword i = 0; i < core_idx.n_elem; ++i) {
      occ_alpha(core_idx(i)) = 1.0;
      occ_beta(core_idx(i))  = 1.0;
    }
  }
  if (nopen > 0) {
    arma::uvec candidates;
    if ((int)order_by_mo.n_elem > ncore) {
      candidates = order_by_mo.subvec(ncore, order_by_mo.n_elem - 1);
    } else {
      candidates.reset();
    }
    if (candidates.n_elem > 0) {
      arma::vec moa_vals(candidates.n_elem);
      for (arma::uword i = 0; i < candidates.n_elem; ++i) moa_vals(i) = moa(candidates(i));
      arma::uvec open_order = arma::sort_index(moa_vals);
      int take = std::min((int)candidates.n_elem, nopen);
      for (int t = 0; t < take; ++t) {
        arma::uword idx = candidates(open_order(t));
        occ_alpha(idx) = 1.0; 
      }
    }
  }
  if ((int)std::floor(spin_n_elec + 0.5) == n_alpha) {
    return occ_alpha;
  } else {
    return occ_beta;
  }
};

scf::DensityMatrix<double> guess = initial_guess_rohf(input, H0, overlap, atoms);
  if (input.at("mixing") != "simple_mixing") {
    throw Error("only simple mixing is implemented in this version");
  }
  const double alpha_mix = input.at("mixing_alpha");
  const auto update_method = scf::simple_mixing<scf::DensityMatrix<double>>{alpha_mix, guess};
  arma::vec n_elec_spin = {(double) n_alpha, (double) n_beta};
  const int max_iter = input.at("max_iter");
  const double energy_tolerance = input.at("energy_tolerance");
  const auto scf_result = scf::scf<scf::simple_mixing<scf::DensityMatrix<double>>, double>(
                          energy_builder,
                          fock_builder,
                          occupation_builder,
                          update_method,
                          overlap,
                          guess,
                          n_elec_spin,
                          max_iter,
                          energy_tolerance);

  nlohmann::json scf_json;
  tr::put(scf_json, "basis_labels", basis.function_label);
  tr::put(scf_json, "energy", scf_result.energy);
  return scf_json;
}





}