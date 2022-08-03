/****************************************************************************
* Copyright (C), 2020-2031 清华大学航天航空学院动力学与控制实验室
* 作者: 张众 zhong-zh19@mails.tsinghua.edu.cn
* 文件名: beam_search.cpp
* 内容简述：通过集束搜索（beam search）实现J2摄动下空间碎片清理任务，源文件
*           后改为地面目标观测目标点问题（CTOC11）
*           代码命名按照Google代码规范
* 文件历史：
* 版本号     日期         作者       说明
* 01        2020-10-03    张众     根据CTOC11修改
* 02        2020-12-17    张众     赛后复盘
****************************************************************************/
#ifndef _PSO_H_
#define _PSO_H_
#include<stdlib.h>
#include<time.h>
#include<iostream>
#include<iomanip>
#include<math.h>
#include<fstream>
#include<random>



//在[min,max)范围内生成随机实数
inline double realRand(const double& min, const double& max) {
	static thread_local std::mt19937_64 generator(std::random_device{}());
	std::uniform_real_distribution<double> distribution(min, max);
	return distribution(generator);
}
//在[min,max]范围内生成随机整数
inline int intRand(const int& min, const int& max) {
	static thread_local std::mt19937_64 generator(std::random_device{}());
	std::uniform_int_distribution<int> distribution(min, max);
	return distribution(generator);
}

/****************************************************************************
* 函数名   : PSO
* 功  能   : PSO优化函数
* 输  入   : X为优化变量，ConstX为参数， xbest为输入量输出为最终结果，fbest为输出指标，D为优化变量个数，Np为种群数，Itmax迭代总数
*            注意：所有变量必须归一化即x的范围是0-1
*            一般种群数为变量个数的10倍，迭代次数1000
****************************************************************************/
void PSO(double (*ObjFun)(double* X, const double* constX), const double* constX, double* xbest, double& fbest, int D, int Np, int ItMax=1000, int ItOut=50,
	double OmegaMin=0.4, double OmegaMax=0.9, double C1Min=0.5, double C1Max=2.5, double C2Min=0.5, double C2Max=2.5, double Vmax=0.3);

void PSO(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* constX ), void* constX, std::vector<double>& xbest, double& fbest, int Np, int ItMax = 1000, int ItOut = 50,
	double OmegaMin = 0.4, double OmegaMax = 0.9, double C1Min = 0.5, double C1Max = 2.5, double C2Min = 0.5, double C2Max = 2.5, double Vmax = 0.3);

#endif