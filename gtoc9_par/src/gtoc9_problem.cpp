/****************************************************************************
* Copyright (C), 2020-2031 Tsinghua University, School of Aerospace Engineering, LAD
* Author: Zhong Zhang
				zhong-zh19@mails.tsinghua.edu.cn
*               539977562@qq.com
* File: gtoc9_problem.cpp
* Description: mission information storage of GTOC9
*
* Log:
*Version      Date        Author           Description
* 01        2022-03-15    Zhong Zhang       Create
****************************************************************************/

#include "gtoc9_problem.h"

#include <iostream>

#include "model.h"

/****************************************************************************
* Function     : mp_calc
* Description  : calculate mp and initial mass based on several velocity increment
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

