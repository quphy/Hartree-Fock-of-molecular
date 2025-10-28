#ifndef hf_rohf_h
#define hf_rohf_h

#include <armadillo>
#include <json.hpp>
#include "../geometry/geometry.hpp"
#include "../basis/basis.hpp"

namespace hf::hf_cpp{
nlohmann::json rohf(const nlohmann::json & input,
                    const geometry::atoms & atoms,
                    const basis::basis & basis);
}
#endif