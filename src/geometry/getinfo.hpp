#ifndef geometry_getinfo
#define geometry_getinfo

#include<string>
#include<armadillo>
#include<vector>
#include<json.hpp>

#include "geometry.hpp"
#include "periodic_table.h"


namespace hf::geometry{
      atoms resolve (const nlohmann::json & inp);
}


#endif