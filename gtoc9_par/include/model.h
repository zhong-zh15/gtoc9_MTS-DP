/****************************************************************************
* Copyright (C), 2020-2031 �廪��ѧ���캽��ѧԺ����ѧ�����ʵ����
* ����: ���� zhong-zh19@mails.tsinghua.edu.cn
*            539977562@qq.com
* �ļ���: model.h
* ���ݼ���������ģ����Ϣ
*
* �ļ���ʷ��
* �汾��     ����         ����       ˵��
* 01a       2021-05-15    ����      �����ļ�
****************************************************************************/
#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <math.h>
#include "Constant.h"
//#include "main.h"

extern double debris_data[123][8];   //ȫ����Ƭ��Ϣ
extern double domega_debris[123];
extern double dOmega_debris[123];
extern double dOmega_debris_notreal[123];

double estimate_T(const std::vector<int>& R, std::vector<double>& T, double initT);
double estimate_dv(int debris_now, int target, double tnow, double& Ts, double& Tf);
double estimate_dv(int debris_now, int target, double tnow, double& Ts, double& Tf, int N);
double estimate_dv(int debris_now, int target, double tnow, double& Ts, double& Tf, double end_epoch, int N);
double Dv_ij(int i, int j, double ts, double tf);
double Dv_All(double* T, int* R, int n);
double Dv_All(const std::vector<double>& T, int* R, int n);
double Dv_All(const std::vector<double>& T, const std::vector<int>& R, std::vector<double>& dv_sequence);
double estimate_dv_random(int debris_now, int target, double tnow, double& Ts, double& Tf, double end_epoch, int N);
double domega(int i);
double domega_init(int i);
double dOmega_init(int i);
double dOmega_notreal(int i);
#endif

