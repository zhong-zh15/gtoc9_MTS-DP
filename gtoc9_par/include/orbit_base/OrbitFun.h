#ifndef _ORBITFUN_H_
#define _ORBITFUN_H_
#include <stdio.h>
#include"Constant.h"

/**************************************************************************************************************************************/
/****************************************************经典轨道根数三种角度关系**********************************************************/
/**************************************************************************************************************************************/
// E2f 根据偏近点角和偏心率求真近点角
double E2f(int& flag, double E, double e);
// E2M 根据偏近点角和偏心率求平近点角
double E2M(int& flag, double E, double e);
// f2E 根据真近点角和偏心率求偏近点角
double f2E(int& flag, double f, double e);
// M2E 根据平近点角和偏心率求偏近点角
double M2E(int& flag, double M, double e, int MaxIter = 100, double epsilon = 1.0e-14);
// f0dt2ft 根据初始真近点角和演化时间求最终真近点角
double f0dt2ft(int& flag, double f0, double dt, double a, double e, double mu, int MaxIter = 100, double epsilon = 1.0e-14);
// f0ft2dt 根据初始真近点角和最终真近点角求演化时间
//double f0ft2dt(int& flag, double f0, double ft, double a, double e, double mu = 3.98600441800e+14, int MaxIter = 100, double epsilon = 1.0e-14);
double f0ft2dt(int& flag, double f0, double ft, double a, double e, double mu = 3.98600441500e+14);
/**************************************************************************************************************************************/
/*******************************************************春分点轨道根数角度关系*********************************************************/
/**************************************************************************************************************************************/
double L0dt2Lt(int& flag, double L0, double dt, const double* ee, double mu);
double L0Lt2dt(int& flag, double L0, double Lt, const double* ee, double mu);
/**************************************************************************************************************************************/
/*********************************************经典轨道根数、直角坐标、改进春分点轨道根数的转换*****************************************/
/**************************************************************************************************************************************/
// coe2rv 根据经典轨道根数求惯性直角坐标系下的位置和速度分量
void coe2rv(int& flag, double* rv, const double* coe, double mu);
// rv2coe 根据地心惯性直角坐标系下的位置和速度分量求经典轨道根数
void rv2coe(int& flag, double* coe, const double* RV, double mu);
// coe2ee 根据经典轨道根数求改进春分点轨道根数，对轨道倾角180度奇异
void coe2ee(int& flag, double* ee, const double* coe, double mu);
// ee2coe 根据改进春分点根数求经典轨道根数
void ee2coe(int& flag, double* coe, const double* ee, double mu);
// ee2rv 根据改进春分点轨道根数求地心惯性直角坐标系下的位置和速度分量
void ee2rv(int& flag, double* rv, const double* ee, double mu);
// rv2ee 根据地心惯性直角坐标系下的位置和速度分量求改进春分点轨道根数，对轨道倾角180度时奇异
void rv2ee(int& flag, double* ee, const double* RV, double mu);

/**************************************************************************************************************************************/
/************************************************已知初值和时间，求末端状态************************************************************/
/**************************************************************************************************************************************/
//根据初始时刻状态coe0求末端时刻dt的状态coe1，按二体推进。若计算成功,flag返回1
void coe02coef(int& flag, double* coe1, const double* coe0, double dt, double mu = 3.98600441500e+14);
//根据初始时刻t0的状态rv0求末端时刻t1的状态rv1，按二体推进。若计算成功,flag返回1
void rv02rvf(int& flag, double* rv1, const double* rv0, double dt, double mu = 3.98600441500e+14);
//根据初始时刻状态ee0求末端时刻dt的状态ee1，按二体推进。若计算成功,flag返回1
void ee02eef(int& flag, double* ee1, const double* ee0, double dt, double mu = 3.98600441500e+14);
//积分求解,根据初始时刻t0的状态rv0求末端时刻t1的状态rv1，按二体推进。若计算成功,flag返回1
void rv02rvf_file(int& flag, double* rv1, const double* rv0, double dt, FILE* fp, double mu = 3.98600441500e+14);




#endif