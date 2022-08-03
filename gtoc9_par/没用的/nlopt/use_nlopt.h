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
* 函数名   : nloptmain()
* 功  能   : nlopt局部优化的调用主函数，优化所有时间
****************************************************************************/
void nloptmain(double (*ObjFun)(const std::vector<double>& x, std::vector<double>& grad, void* f_data), int* node_order,
	std::vector<double>& xbest, double& fbest, int ItMax);

void map_input(vector<double>& X);
/****************************************************************************
* 函数名   : nlopt_tstf
* 功  能   : 优化的目标函数；返回值是优化指标
****************************************************************************/
double nlopt_tstf(const std::vector<double>& X, std::vector<double>& grad, void* f_data);
double nlopt_tstf_pso(const std::vector<double>& x, std::vector<double>& grad, void* f_data);
#endif
