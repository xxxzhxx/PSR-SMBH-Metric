#include "metric.h"

#include <math.h>
#include <stdio.h>

void metric_guu_Schwarzschild(const double xu[4], double guu[4][4], void *params){
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            guu[i][j]=0;
        }
    }

    double M=*(double *)params;
    double r=xu[1],theta=xu[2];
    double sin_theta=sin(theta);
    double C=M/r;
    double g=1.0-2.0*C;

    guu[0][0]=-1.0/g;
    guu[1][1]=g;
    guu[2][2]=1.0/(r*r);
    guu[3][3]=1.0/(r*r*sin_theta*sin_theta);
}
void metric_gdd_Schwarzschild(const double xu[4], double gdd[4][4], void *params){
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            gdd[i][j]=0;
        }
    }

    double M=*(double *)params;
    double r=xu[1],theta=xu[2];
    double sin_theta=sin(theta);
    double C=M/r;
    double g=1.0-2.0*C;

    gdd[0][0]=-g;
    gdd[1][1]=1.0/g;
    gdd[2][2]=r*r;
    gdd[3][3]=r*r*sin_theta*sin_theta;
}
void affine_Gudd_Schwarzschild(const double xu[4], double Gudd[4][4][4], void *params){
    double M=*(double *)params;
    double r=xu[1],theta=xu[2];
    double sin_theta=sin(theta),cos_theta=cos(theta);
    double C=M/r;
    double g=1.0-2.0*C;

    Gudd[0][0][0]=0;
    Gudd[0][0][1]=C/r/g;
    Gudd[0][0][2]=0;
    Gudd[0][0][3]=0;
    Gudd[0][1][1]=0;
    Gudd[0][1][2]=0;
    Gudd[0][1][3]=0;
    Gudd[0][2][2]=0;
    Gudd[0][2][3]=0;
    Gudd[0][3][3]=0;
    
    Gudd[1][0][0]=C*g/r;
    Gudd[1][0][1]=0;
    Gudd[1][0][2]=0;
    Gudd[1][0][3]=0;
    Gudd[1][1][1]=0;
    Gudd[1][1][2]=-C/r/g;
    Gudd[1][1][3]=0;
    Gudd[1][2][2]=-r*g;
    Gudd[1][2][3]=0;
    Gudd[1][3][3]=-r*g*sin_theta*sin_theta;

    Gudd[2][0][0]=0;
    Gudd[2][0][1]=0;
    Gudd[2][0][2]=0;
    Gudd[2][0][3]=0;
    Gudd[2][1][1]=0;
    Gudd[2][1][2]=1.0/r;
    Gudd[2][1][3]=0;
    Gudd[2][2][2]=0;
    Gudd[2][2][3]=0;
    Gudd[2][3][3]=-sin_theta*cos_theta;

    Gudd[3][0][0]=0;
    Gudd[3][0][1]=0;
    Gudd[3][0][2]=0;
    Gudd[3][0][3]=0;
    Gudd[3][1][1]=0;
    Gudd[3][1][2]=0;
    Gudd[3][1][3]=1.0/r;
    Gudd[3][2][2]=0;
    Gudd[3][2][3]=cos_theta/sin_theta;
    Gudd[3][3][3]=0;

    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            for(int k=0;k<j;k++){
                Gudd[i][j][k]=Gudd[i][k][j];
            }
        }
    }
}

