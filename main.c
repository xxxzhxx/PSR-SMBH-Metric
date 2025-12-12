#include <stdio.h>
#include <time.h>
#include <math.h>

#include <gsl/gsl_math.h>

#include "metric.h"
#include "motion.h"

int main(int argc, char *argv[]){
    double total_t=0;
    clock_t start_t,finish_t;
    start_t=clock();

    struct metric m;
    m.guu=metric_guu_Schwarzschild;
    m.gdd=metric_gdd_Schwarzschild;
    m.Gudd=affine_Gudd_Schwarzschild;
    double M=1.0;

    struct motion_system mot_sys;
    mot_sys.m=m;
    mot_sys.metric_params=&M;
    mot_sys.massive=1;
    mot_sys.tau0=0.0;

    double r=100.0*M;
    mot_sys.rv0[0]=0.0;
    mot_sys.rv0[1]=r;
    mot_sys.rv0[2]=M_PI_2;
    mot_sys.rv0[3]=0.0;
    mot_sys.rv0[5]=0.0;
    mot_sys.rv0[6]=0.0;
    mot_sys.rv0[7]=sqrt(M/r)/sqrt(1.0-3.0*M/r)/r*1.0/2.0;
    mot_sys.rv0[4]=sqrt(1.0+r*r*mot_sys.rv0[7]*mot_sys.rv0[7])/sqrt(1.0-2.0*M/r);

    int n=3000;
    double tau[n];
    for(int i=0;i<n;i++){
        tau[i]=(double)i*1.0;
    }
    struct motion_output out[n];
    motion_evolve_constraint(mot_sys,n,tau,out);
    FILE *fp=fopen("check_code/motion_constraint.txt","wt");
    for(int i=0;i<n;i++){
        fprintf(fp,"%.16g ",out[i].tau);
        for(int j=0;j<8;j++){
            fprintf(fp,"%.16g ",out[i].rv[j]);
        }
        fprintf(fp,"%.16g\n",out[i].error_track);
    }
    fclose(fp);

    
    finish_t=clock();
    total_t=(double)(finish_t-start_t)/CLOCKS_PER_SEC;
    printf("num calculation in %g sec\n",total_t);
    return 0;
}