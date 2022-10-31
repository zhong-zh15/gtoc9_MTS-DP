/****************************************************************************
* Copyright (C), 2020-2031 Tsinghua University, School of Aerospace Engineering, LAD
* Author: Zhong Zhang
				zhong-zh19@mails.tsinghua.edu.cn
*               539977562@qq.com
* File: gtoc9_problem.h
* Description: mission information storage of GTOC9
*
* Log:
*Version      Date        Author           Description
* 01        2022-03-15    Zhong Zhang       Create
****************************************************************************/
#pragma once
#include <vector>

#include "Constant.h"

#include "DP.h"

using namespace std;

/****************************************************************************
* Function     : mp_calc
* Description  : calculate mp and initial mass based on several velocity increment
****************************************************************************/
double mp_calc(const std::vector<double>& dv, double& m0);

/****************************************************************************
* Struct       : Opt_info_gtoc9
* Description  : used to save sequence info
****************************************************************************/
struct Opt_info_gtoc9
{
	int mission;
	double end_peoch;
	int last_debris_num;
	vector<vector<double>> time_sequence;
	vector<vector<int>> debris_squence;
};

/****************************************************************************
* Struct       : opt_struct
* Description  : save each mission's info
****************************************************************************/
struct opt_struct
{
	double  optimin;
	double opt_each_mission[TreeNum];
	std::vector<vector<double>> t;
	std::vector<vector<double>> dv;
};