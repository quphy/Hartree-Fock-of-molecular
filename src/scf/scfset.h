#ifndef setup_h
#define setup_h

#include "scf.h"
#include "../geometry/geometry.hpp"
#include "../basis/basis.hpp"

namespace hf::scf{
template<typename T>
struct Setup{
  geometry::atoms atoms;
  basis::basis basis;
  scf::OverlapMatrix<T> overlap;
};

template<class Setup,typename T>
struct Result{
  Setup setup;
  scfresult<T> scf_result;
};
}


#endif