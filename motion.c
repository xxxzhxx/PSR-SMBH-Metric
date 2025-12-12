#include "motion.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include "metric.h"

static int func_motion_evolve(double tau, const double rv[], double f[], void *params){
    //rv里是：x0,x1,x2,x3,dx0,dx1,dx2,dx3
    struct motion_system *m=(struct motion_system *)params;
    double xu[4]={rv[0],rv[1],rv[2],rv[3]};
    double Gudd[4][4][4];
    m->m.Gudd(xu,Gudd,m->metric_params);

    // 考虑每步给v归一化一下？
    // 但有个问题是，对于类光和类时是不一样的，类光没法简单说归一化。。
    // 硬降一维？
    // 或者。。。加一维？让error_track进入到积分步长控制里去
    // 这可能需要使用更底层的控制

    f[0]=rv[4];
    f[1]=rv[5];
    f[2]=rv[6];
    f[3]=rv[7];

    f[4]=0;
    f[5]=0;
    f[6]=0;
    f[7]=0;

    for(int i=0;i<4;i++){
        f[4]+=-Gudd[0][i][i]*f[i]*f[i];
        f[5]+=-Gudd[1][i][i]*f[i]*f[i];
        f[6]+=-Gudd[2][i][i]*f[i]*f[i];
        f[7]+=-Gudd[3][i][i]*f[i]*f[i];
    }

    for(int i=0;i<4;i++){
        for(int j=i+1;j<4;j++){
            f[4]+=-2.0*Gudd[0][i][j]*f[i]*f[j];
            f[5]+=-2.0*Gudd[1][i][j]*f[i]*f[j];
            f[6]+=-2.0*Gudd[2][i][j]*f[i]*f[j];
            f[7]+=-2.0*Gudd[3][i][j]*f[i]*f[j];
        }
    }

    return GSL_SUCCESS;
}
//注意，我们强制要求系统处于tau[0]时刻，否则产生理解的歧义（由于tau只是参数，其实不影响），后果自负
//这里明显第一个结构体要改一下
void motion_evolve(struct motion_system mot_sys, int n, double tau[], struct motion_output out[]){
    if(mot_sys.tau0!=tau[0]){
        //直接给个报错抬走
        printf("motion_evolve error: please check the inital tau0 and tau[0]\n");
        return;
    }

    gsl_odeiv2_system sys={func_motion_evolve,NULL,8,&mot_sys};
    gsl_odeiv2_driver *d=gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rk8pd,1e-8,1e-14,1e-14);

    double tau_now=tau[0];
    double rv[8];
    for(int i=0;i<8;i++){
        rv[i]=mot_sys.rv0[i];
    }

    for(int i=0;i<n;i++){
        gsl_odeiv2_driver_apply(d,&tau_now,tau[i],rv);
        out[i].tau=tau_now;
        for(int j=0;j<8;j++){
            out[i].rv[j]=rv[j];
        }
        double xu[4]={rv[0],rv[1],rv[2],rv[3]};
        double gdd[4][4];
        mot_sys.m.gdd(xu,gdd,mot_sys.metric_params);
        out[i].error_track=0;
        for(int j=0;j<4;j++){
            for(int k=0;k<4;k++){
                out[i].error_track+=gdd[j][k]*rv[4+j]*rv[4+k];
            }
        }
        if(mot_sys.massive==1){
            out[i].error_track+=1.0;
        }
    }

    gsl_odeiv2_driver_free(d);
    return;
}

void *my_control_alloc(void){
    //不需要额外空间
    return NULL;
}
static void my_control_free(void *state){
    //不需要释放mot_sys
    return;
}
static int my_control_init(void *state, double eps_abs, double eps_rel, double a_y, double a_dydt){
    //没有初始化
    return GSL_SUCCESS;
}
static int my_control_hadjust(void *state, size_t dim, unsigned int ord, const double y[], const double yerr[], const double yp[], double *h){
    //这是核心，我们要实现的功能是处理以下三种情况
    /* step was increased */
    /* step unchanged     */
    /* step decreased     */
    struct motion_system *mot_sys=(struct motion_system *)state;
    double violation=0;
    double xu[4]={y[0],y[1],y[2],y[3]};
    double gdd[4][4];
    mot_sys->m.gdd(xu,gdd,mot_sys->metric_params);
    for(int i=0;i<4;i++){
        violation+=gdd[i][i]*y[4+i]*y[4+i];
    }
    for(int i=0;i<4;i++){
        for(int j=0;j<i;j++){
            violation+=2.0*gdd[j][i]*y[4+j]*y[4+i];
        }
    }
    if(mot_sys->massive==1){
        violation+=1.0;
    }

    double err_tol=1e-4;
    if(violation>err_tol){
        //应该缩步长
        *h *=0.9;
        return GSL_ODEIV_HADJ_DEC;
    }else if(violation<0.5*err_tol){
        //可以放步长
        *h /=0.9;
        printf("%.16g %.16g %.16g\n",y[0],violation,*h);
        return GSL_ODEIV_HADJ_INC;
    }else{
        //接受结果并保持步长不变
        printf("%.16g %.16g %.16g\n",y[0],violation,*h);
        return GSL_ODEIV_HADJ_NIL;
    }
}
int my_control_errlevel(void *state, const double y, const double dydt, const double h, const size_t ind, double *errlev){
    //不是传统控制
    *errlev=1.0;
    return GSL_SUCCESS;
}
int my_control_set_driver(void *state, const gsl_odeiv2_driver * d){
    //这个函数据说是用driver来调节state，因此我们不需要
    return GSL_SUCCESS;
}
void motion_evolve_constraint(struct motion_system mot_sys, int n, double tau[], struct motion_output out[]){
    if(mot_sys.tau0!=tau[0]){
        //直接给个报错抬走
        printf("motion_evolve error: please check the inital tau0 and tau[0]\n");
        return;
    }

    gsl_odeiv2_system sys={func_motion_evolve,NULL,8,&mot_sys};
    gsl_odeiv2_driver *d=gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rk8pd,1e-8,1e-14,1e-14);

    //我们用自己定义的步长控制取代内置控制
    gsl_odeiv2_control_free(d->c);
    gsl_odeiv2_control c;
    d->c=&c;

    c.state=&mot_sys;
    gsl_odeiv2_control_type control_type;
    c.type=&control_type;

    control_type.name="my_control";
    control_type.alloc=my_control_alloc;
    control_type.init=my_control_init;
    control_type.hadjust=my_control_hadjust;
    control_type.errlevel=my_control_errlevel;
    control_type.set_driver=my_control_set_driver;
    control_type.free=my_control_free;

    double tau_now=tau[0];
    double rv[8];
    for(int i=0;i<8;i++){
        rv[i]=mot_sys.rv0[i];
    }

    for(int i=0;i<n;i++){
        gsl_odeiv2_driver_apply(d,&tau_now,tau[i],rv);
        out[i].tau=tau_now;
        for(int j=0;j<8;j++){
            out[i].rv[j]=rv[j];
        }
        double xu[4]={rv[0],rv[1],rv[2],rv[3]};
        double gdd[4][4];
        mot_sys.m.gdd(xu,gdd,mot_sys.metric_params);
        out[i].error_track=0;
        for(int j=0;j<4;j++){
            for(int k=0;k<4;k++){
                out[i].error_track+=gdd[j][k]*rv[4+j]*rv[4+k];
            }
        }
        if(mot_sys.massive==1){
            out[i].error_track+=1.0;
        }
    }

    gsl_odeiv2_driver_free(d);
    return;
}
