#ifndef print_h
#define print_h
#include <armadillo>
#include <fmt/format.h>



namespace hf{


//针对数字转化为字符串string的模版
template <typename T>
std::string format(const T number,const int precision=6,
                   const int width=-1,const std::string aligned=">" ){

if(width<0){
    return fmt::format("{}",number);
}
else {
    //{:.{}}指定小数点后保留的位数
    const auto formatted=fmt::format("{:.{}}",number,precision);
    return fmt::format("{:"+aligned+"{}}",formatted,width);   
}
}

//针对double的特化
template<>
inline
std::string format(const double number,const int precision,
                   const int width,const std::string aligned) {
  const auto formatted=fmt::format("{:.{}g}",number,precision);
  if(width < 0) {
    return formatted;
  } 
  else {
    return fmt::format("{:"+aligned+"{}}",formatted,width);
  }
}

//针对复数的特化
template<>
inline
std::string format(const arma::cx_double number,const int precision,
                   const int width,const std::string aligned) {
  const auto formatted ="("+fmt::format("{:.{}g}",std::real(number),precision)+","+
                         fmt::format("{:.{}g}",std::imag(number),precision)+")";

  if(width < 0) {
    return formatted;
  } 
  else {
    return fmt::format("{:"+aligned +"{}}",formatted,width);
  }
}

//将转化的字符串输出的模版
template <typename T>
void print(const arma::Row<T>&arma,const int width =10,
           const int precison=6, const std::string aligned=">" ){
  for (arma::uword j=0;j<arma.n_cols;j++){
    const std::string formatted=format<T>(arma(j),precison,width,aligned);
    fmt::print("{:"+aligned+"{}}",formatted,width);
  }
}

template <typename T>
void print(const arma::Col<T>&arma,const int width =10,
           const int precison=6, const std::string aligned=">" ){
  for (arma::uword j=0;j<arma.n_rows;j++){
    const std::string formatted=format<T>(arma(j),precison,width,aligned);
    fmt::print("{:"+aligned+"{}}",formatted,width);
    fmt::print("\n");
  }
}

template <typename T>
void print(const arma::Mat<T>&arma, const int width=10,
           const int precision=6,const std::string aligned=">"){
  for (arma::uword i=0; i<arma.n_rows;i++) {
    for (arma::uword j=0; j<arma.n_cols;j++) {
      const std::string formatted =format<T>(arma(i,j), precision,width,aligned);
      fmt::print("{:"+aligned+"{}}",formatted,width);
    }
    fmt::print("\n");
  }
}

template <typename T>
void print (const T number,const int width=10,
            const int precision =6,const std::string aligned=">"){
  const auto formatted = format<T>(number, precision,width,aligned);
  fmt::print("{:" + aligned + "{}}", formatted, width);
}


template<typename State>
using Printer = std::function<int(const State & state,
                                  const double time,
                                  const int iter)>;

//========分裂符号========
inline
void print_separator(const int length, const std::string symbol = "=") {
  std::string line_separator = "";
  for (int i = 0; i < length; i++) {
    line_separator += symbol;
  }
  fmt::print(line_separator + "\n");
}
}

#endif