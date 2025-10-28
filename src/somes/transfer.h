#ifndef transfer_json_armadillo
#define transfer_json_armadillo

#include <iostream>
#include <vector>
#include <armadillo>
#include <json.hpp>

#include "error&warning.hpp"

//使用nlohamnn的json，参考官方文档https://json.nlohmann.me/api/basic_json/

namespace hf {
   namespace tr{

  //template function NO.2
  //json转化成数组
    template<typename T>
 //template函数，将一个json数组转变为cpp向量
    std::vector<T> get_vecor_from_json(const nlohmann::json & jsonarray) {
      std::vector<T> result;
 //遍历数组元素，get是json函数
      for (const auto & jsonelement : jsonarray) {
        result.push_back(jsonelement.get<T>());
      }
      return result;
    }
    
    //template function NO.3 
    //提取完行列互换，但是符合armadilo和cpp列优先
    template<typename T>
    arma::Mat<T> get_mat_from_json(const nlohmann::json & jsonmatric) {
      arma::Mat<T> result;
    //遍历matric每一行
      for (const auto & jsonline : jsonmatric) {
        //提取每一行，转置成列
        const auto jsonelements = arma::Col<T>(get_vecor_from_json<T>(jsonline));
        //拼成armadillo矩阵
        result = arma::join_rows(result, jsonelements);
      }
      return result;
    }
    

  
   //template function NO.1 
   //模板函数，T为抽象概念 
     template<typename T>
    T get(const nlohmann::json & jsonabstract) {
      return jsonabstract.get<T>();
    }

    //特化NO.1函数，针对arma::vec，arma::cx_vec，arma::mat，arma::cx_mat，arma::cube数据类型，
    //arma::vec列向量，arma::cx_vec复数列向量，arma::mat矩阵，arma::cx_mat复数矩阵
    //arma::cube三维数组（立方体）
    template<>
    inline
    arma::vec get(const nlohmann::json & jsonabstract) {
      return get_vecor_from_json<double>(jsonabstract);
    }
    
    template<>
    inline
    arma::cx_vec get(const nlohmann::json & jsonabstract) {
      const arma::vec real = get_vecor_from_json<double>(jsonabstract["real"]);
      const arma::vec imag = get_vecor_from_json<double>(jsonabstract["imag"]);
      return arma::cx_vec{real, imag};
    }
    
    template<>
    inline
    arma::mat get(const nlohmann::json & jsonabstract) {
      return get_mat_from_json<double>(jsonabstract);
    }
    
    template<>
    inline
    arma::cx_mat get(const nlohmann::json & jsonabstract) {
      const arma::mat real = get_mat_from_json<double>(jsonabstract["real"]);
      const arma::mat imag = get_mat_from_json<double>(jsonabstract["imag"]);
      return arma::cx_mat{real, imag};
    }
    
    template<>
    inline
    arma::cube get(const nlohmann::json & jsonabstract) {
      const auto length = jsonabstract.size();
      const arma::mat sample = get<arma::mat>(jsonabstract[0]);
      arma::cube result(sample.n_rows, sample.n_cols, length);
      for (unsigned long i=0; i<length; i++) {
        result.slice(i) = get<arma::mat>(jsonabstract[i]);
      }
      return result;
    }
    
    //template function no4
    //数组转化成json
    //模板函数
    template<typename T>
    auto put(nlohmann::json & result, const std::string & path, const T & value) {
      (result[path] = value);
    //put(result, "data", 25)的结果是，result中储存了
    /*
    {
    "data": 25
    }
    */
    }
    //特化
    template<>
    inline
    auto
    put(nlohmann::json & result, const std::string & path, const arma::vec & value) {
    
      std::vector<double> converted_array;
    
      for(arma::uword i=0; i<value.n_elem; i++) {
        converted_array.push_back(value(i));
      }
    
      result[path] = converted_array;
    }
    
    template<>
    inline
    auto
    put(nlohmann::json & result, const std::string & path, const arma::mat & value) {
    
      std::vector<std::vector<double>> converted_matrix;
    
      for(arma::uword i=0; i<value.n_rows; i++) {
        std::vector<double> row_vector;
        for(arma::uword j=0; j<value.n_cols; j++) {
          row_vector.push_back(value(i, j));
        }
        converted_matrix.push_back(row_vector);
      }
    
      (result[path] = converted_matrix);
    }
    
    template<>
    inline
    auto put(nlohmann::json & result, const std::string & path,
             const arma::cx_vec & value) {
    
      nlohmann::json head;
      const arma::vec real = arma::real(value);
      const arma::vec imag = arma::imag(value);
      put(head, "real", real);
      put(head, "imag", imag);
    
      result[path] = head;
    }
    
    template<>
    inline
    auto put(nlohmann::json & result, const std::string & path,
             const arma::cx_mat & value) {
    
      nlohmann::json head;
    
      const arma::vec real = arma::real(value);
      const arma::vec imag = arma::imag(value);
    
      put(head, "real", real);
      put(head, "imag", imag);
    
      result[path] = head;
    }
    
    
    template<typename T>
    inline
    auto put(nlohmann::json & result, const std::string & path,
             const std::vector<T> & value) {
      result[path] = value;
    }
    
    }
    }

#endif