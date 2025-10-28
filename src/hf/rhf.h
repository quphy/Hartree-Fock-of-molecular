#ifndef hf_rhf_h
#define hf_rhf_h

#include <json.hpp>

#include "../geometry/geometry.hpp"
#include "../basis/basis.hpp"
#include "../scf/scf.h"
#include "../scf/scfset.h"

namespace hf::hf_cpp{
template<typename T>
struct RHFSetup : scf::Setup<T> {
  arma::mat H0;
  scf::FockBuilder<T> fock_builder;
  scf::EnergyBuilder<T> energy_builder;
  scf::OccupationBuilder occupation_builder;
  scf::DensityMatrix<T> initial_guess;
};

template<typename T>
using RHFResult=scf::Result<RHFSetup<T>, T>;

arma::mat core_hamiltonian(const geometry::atoms & atoms, const basis::basis & basis);

nlohmann::json rhf(const nlohmann::json & input,const geometry::atoms & atoms,const basis::basis & basis);

}
#endif