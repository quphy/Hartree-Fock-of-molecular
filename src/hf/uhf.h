#ifndef hf_uhf_h
#define hf_uhf_h
#include <json.hpp>
#include "geometry/geometry.hpp"
#include "basis/basis.hpp"
#include "scf/scf.h"
#include "scf/scfset.h"

namespace hf::hf_cpp{

template<typename T>
struct UHFSetup : scf::Setup<T> {
  arma::mat H0;
  scf::FockBuilder<T> fock_builder;
  scf::EnergyBuilder<T> energy_builder;
  scf::OccupationBuilder occupation_builder;
  scf::DensityMatrix<T> initial_guess;
};

template<typename T>
using UHFResult = scf::Result<UHFSetup<T>, T>;

nlohmann::json uhf(const nlohmann::json &input,const geometry::atoms & atoms,const basis::basis & basis);


}






#endif