#ifndef  USENLOPT
#define  USENLOPT


#include <vector>
#include <iterator>

#include "nlopt.hpp"
#include "model.h"

using std::vector;
using std::begin;
using std::end;
/****************************************************************************
* ������   : nloptmain()
* ��  ��   : nlopt�ֲ��Ż��ĵ������������Ż�����ʱ��
****************************************************************************/
void nloptmain(double (*ObjFun)(const std::vector<double>& x, std::vector<double>& grad, void* f_data), int* node_order,
	std::vector<double>& xbest, double& fbest, int ItMax);

void map_input(vector<double>& X);
/****************************************************************************
* ������   : nlopt_tstf
* ��  ��   : �Ż���Ŀ�꺯��������ֵ���Ż�ָ��
****************************************************************************/
double nlopt_tstf(const std::vector<double>& X, std::vector<double>& grad, void* f_data);
double nlopt_tstf_pso(const std::vector<double>& x, std::vector<double>& grad, void* f_data);
#endif
