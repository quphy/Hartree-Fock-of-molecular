#include<string>
#include<armadillo>
#include<vector>
#include<json.hpp>

#include "geometry.hpp"
#include "periodic_table.h"


namespace hf::geometry {
    //将input文件中geometry块中的信息解析到结构体atoms中
    atoms resolve (const nlohmann::json & inp){
        std::vector<std::string> symbols;
        std::vector<arma::uword> atom_numbers;
        std::vector<double> coordinates;
        //电荷读取
        const int charge = inp.at("charge");
        //原子信息读取，检查原子信息，应该有4行
        const auto atom_info = inp.at("atoms");
        for (const auto & line : atom_info) {
            int unit_index = 0;

        

        //通过unit_index，将信息获取，json示例为["H",0,0,0]
            for(const auto & unit : line) {
              if(unit_index == 0) {
                const auto atom_symbol= unit.get<std::string>();
                symbols.push_back(atom_symbol);
                atom_numbers.push_back(geometry::periodic_table.at(atom_symbol));
              } else {
                //读取对应原子的坐标
                const auto coordinate = unit.get<double>();
                coordinates.push_back(coordinate);
              }
              unit_index++;
            }
          }
        //将angstrom转化为bohr，并将长vector 转化为3个一组的坐标，3*symbols.size
          const double bohr = 1.8897259886;
          const arma::mat xyz =
              arma::reshape(arma::vec(coordinates), 3, symbols.size())
              * bohr;
        
          return {symbols, arma::uvec{atom_numbers}, xyz, charge};
 }
}



