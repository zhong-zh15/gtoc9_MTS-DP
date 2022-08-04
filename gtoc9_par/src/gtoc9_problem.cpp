#include "gtoc9_problem.h"

#include <iostream>


#include "model.h"



/****************************************************************************
* 函数名   : mp_calc()
* 功  能   : 根据几次速度增量 计算mp，和初始质量
****************************************************************************/
double mp_calc(const std::vector<double>& dv, double& m0)
{
	double mf = 2000; //kg
	double m_temp = mf;
	double mp = 0.0;
	for (int i = dv.size() - 1; i > -1; i--)
	{
		m_temp += 30;
		double mp_temp = m_temp * exp(dv[i] / 340.0 / 9.80665) - m_temp;
		m_temp = m_temp + mp_temp;
		mp += mp_temp;
	}
	m_temp += 30;
	m0 = m_temp;
	return mp;
}

