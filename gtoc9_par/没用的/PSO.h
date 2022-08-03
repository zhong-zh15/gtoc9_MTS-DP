/****************************************************************************
* Copyright (C), 2020-2031 �廪��ѧ���캽��ѧԺ����ѧ�����ʵ����
* ����: ���� zhong-zh19@mails.tsinghua.edu.cn
* �ļ���: beam_search.cpp
* ���ݼ�����ͨ������������beam search��ʵ��J2�㶯�¿ռ���Ƭ��������Դ�ļ�
*           ���Ϊ����Ŀ��۲�Ŀ������⣨CTOC11��
*           ������������Google����淶
* �ļ���ʷ��
* �汾��     ����         ����       ˵��
* 01        2020-10-03    ����     ����CTOC11�޸�
* 02        2020-12-17    ����     ������
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



//��[min,max)��Χ���������ʵ��
inline double realRand(const double& min, const double& max) {
	static thread_local std::mt19937_64 generator(std::random_device{}());
	std::uniform_real_distribution<double> distribution(min, max);
	return distribution(generator);
}
//��[min,max]��Χ�������������
inline int intRand(const int& min, const int& max) {
	static thread_local std::mt19937_64 generator(std::random_device{}());
	std::uniform_int_distribution<int> distribution(min, max);
	return distribution(generator);
}

/****************************************************************************
* ������   : PSO
* ��  ��   : PSO�Ż�����
* ��  ��   : XΪ�Ż�������ConstXΪ������ xbestΪ���������Ϊ���ս����fbestΪ���ָ�꣬DΪ�Ż�����������NpΪ��Ⱥ����Itmax��������
*            ע�⣺���б��������һ����x�ķ�Χ��0-1
*            һ����Ⱥ��Ϊ����������10������������1000
****************************************************************************/
void PSO(double (*ObjFun)(double* X, const double* constX), const double* constX, double* xbest, double& fbest, int D, int Np, int ItMax=1000, int ItOut=50,
	double OmegaMin=0.4, double OmegaMax=0.9, double C1Min=0.5, double C1Max=2.5, double C2Min=0.5, double C2Max=2.5, double Vmax=0.3);

void PSO(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* constX ), void* constX, std::vector<double>& xbest, double& fbest, int Np, int ItMax = 1000, int ItOut = 50,
	double OmegaMin = 0.4, double OmegaMax = 0.9, double C1Min = 0.5, double C1Max = 2.5, double C2Min = 0.5, double C2Max = 2.5, double Vmax = 0.3);

#endif