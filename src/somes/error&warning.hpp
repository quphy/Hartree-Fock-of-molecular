#ifndef somes_error_h
#define somes_error_h

#include <string>
#include <stdexcept>
#include <iostream>


namespace hf {
//继承自标准库的error类
struct Error : public std::runtime_error {
//构造函数，接受string型变量
  explicit Error(const std::string & error_message) :
      std::runtime_error("Error: " + error_message) {}
//构造函数，接受char 字符串型变量
  explicit Error(const char * error_message) :
      std::runtime_error(std::string("Error: ")+ std::string(error_message)) {}
};

}
#endif 