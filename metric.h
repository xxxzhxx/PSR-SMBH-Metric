#ifndef METRIC_H
#define METRIC_H

// 我们在这里应当定义一些函数，metric 和 affine，示例如下
// 我们强制要求是4d时空，第一维是时间
// void *metric_guu(const double xu[4], double guu[4][4], void *params); 方便起见g不妨成对定义，当然最后用到哪个定义哪个也行
// void *metric_gdd(const double xu[4], double gdd[4][4], void *params);
// void *affine_Gudd(const double xu[4], double Gudd[4][4][4], void *params);
// 实际使用的话怎么方便怎么来，这只是说给一个通用接口，肯定不能做到优化和特化的

//我们传参数的时候应该用结构体，结构体里应该是函数指针。。。
struct metric{
    void (*guu)(const double[4],double[4][4],void *);
    void (*gdd)(const double[4],double[4][4],void *);
    void (*Gudd)(const double[4],double[4][4][4],void *);
};

// 史瓦西度规
void metric_guu_Schwarzschild(const double xu[4], double guu[4][4], void *params);
void metric_gdd_Schwarzschild(const double xu[4], double gdd[4][4], void *params);
void affine_Gudd_Schwarzschild(const double xu[4], double Gudd[4][4][4], void *params);



#endif