#include<iostream>
#include<armadillo>
#include<args.hxx>
#include<fmt/core.h>
#include<json.hpp>
#include<fstream>

#include "somes/time.h"
#include "somes/error&warning.hpp"
#include "cal.hpp"

/*
This project was developed with reference to and by adapting the open-source work of Rui Li, 
originally released in 2020 under the MIT License, located in https://github.com/Walter-Feng/Hartree-Fock-in-CPP.git.
In the process of studying, borrowing and rewriting his code, I learned good C++ Objected Oriented Programming 
and gained a deeper understanding of the specific implementation of Hartree-Fock. 
I hereby express my sincere gratitude to Rui Li for his generous sharing. His project provided essential code resources 
and design ideas that significantly supported the training and implementation of this project.
*/

int main(int argc,char* argv[]){
using namespace hf;
timecounting global_time;
std::cerr.rdbuf(std::cout.rdbuf());
//记录程序总耗时 
//命令行参数解析
//说明
args::ArgumentParser annoum("This is a quantum chemistry program written by Ziqi Wang. To run the program"
                            "./{executable filename} {input file path} >{output file path} "   
                           );  
//help选项
args::HelpFlag help(annoum, "help","Display this help menu", {'h', "help"});
//input格式解释
args::Positional<std::string> inputfile_flag (annoum,"inputfile",
                                        "the required inputfile should be writen in json format");                                            
 
//利用try catch语句分辨./是否有问题
try {
    annoum.ParseCLI(argc, argv);
}
catch (const args::Help &) {
    std::cout << annoum << std::endl;
    return 0;
}
//数据解析失败提示
catch (const args::ParseError & e) {
    std::cout << e.what() << std::endl;
    std::cout << annoum << std::endl;
    return 1;
}
//数据提交失败提示
catch (const args::ValidationError & e) {
    std::cout << e.what() << std::endl;
    std::cout << annoum << std::endl;
    return 1;
}

//读取input文件
std::cout<<"It seems there's nothing wrong, Let's start reading the input"<<std::endl;
nlohmann::json inputfile;
//获取运行的时候传入的文件路径，args::get(inputfile_flag)得到的就是输入文件的路径
const std::string inputfile_filename=args::get(inputfile_flag);
//如果没有给出输入文件
 if(inputfile_filename.empty()) {
    throw Error("File input is not given, use -h to see the help.");
  }
  std::ifstream input_file_stream(args::get(inputfile_flag));

  input_file_stream >> inputfile;
  input_file_stream.close();
 

//利用cal函数进行计算
nlohmann::json result=cal(inputfile);


std::cout<<"HF calculation is done. wish you feel happy"<<std::endl;
    return 0;
}