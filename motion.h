#ifndef MOTION_H
#define MOTION_H

#include <stdbool.h>

#include "metric.h"

struct motion_system{
    struct metric m;
    void *metric_params;
    double tau0;
    double rv0[8];
    bool massive;
};

struct motion_output{
    double tau;
    double rv[8];
    double error_track;
};

// 我们尝试了约束方案，但失败了，所以这里还是采取消去一维，考虑到符号问题，消去dt/dtau是最自然的（其实在轴对称问题里消去dphi/dtau也是自然的？）
void motion_evolve(struct motion_system mot_sys, int n, double tau[], struct motion_output out[]);
// 我们由归一化条件写一个步长控制，先当初步尝试，未必好用，类时类光也未必通用
// 确实不好用，这个东西对步长并不敏感，所以会直接步长掉到很小难以完成积分。。
void motion_evolve_constraint(struct motion_system mot_sys, int n, double tau[], struct motion_output out[]);

#endif