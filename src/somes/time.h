#ifndef some_time_h
#define some_time_h

//用chrono的库定义一个计时器
#include<chrono>

namespace hf{
//定义一个类进行记时，用的时候 time “对象名”，“对象名”.function
class timecounting
{
//利用private保护时钟和单位，以防被篡改
private:
//定义一个最高精度时钟
using defclock = std::chrono::high_resolution_clock;
//利用模板类定义时间间隔以秒和分钟单位
using defsecond = std::chrono::duration<double,std::ratio<1>>;
using defminute = std::chrono::duration<double,std::ratio<60>>;

//用定义的高精度时钟上的时间表示当前时刻
//定义时间点对象start_counting_time，用来储存time类创造的时间点
std::chrono::time_point<defclock> start_counting_time;


//利用public调用private的内容，方便记时
public:
//构造函数，当别处构造time类的时候，那个时刻的就被记录在now_time
timecounting()
{
    start_counting_time = defclock::now();
}

//返回现在的时间，分别以秒和分钟计算
double timesec() const
{
double sec;
sec=std::chrono::duration_cast<defsecond>(defclock::now() - start_counting_time).count();
return sec;
}

double timemin() const
{
double min;
min=std::chrono::duration_cast<defminute>(defclock::now() - start_counting_time).count();
return min;
}

};

}

#endif