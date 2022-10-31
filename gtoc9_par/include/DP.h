/****************************************************************************
* Copyright (C), 2020-2031 Tsinghua University, School of Aerospace Engineering, LAD
* Author: Zhong Zhang
				zhong-zh19@mails.tsinghua.edu.cn
*               539977562@qq.com
* File: DP.h
* Description: the specific implementation of dynamic programming on the GTOC9 problem
*
* Log:
*Version      Date        Author           Description
* 01        2022-03-15    Zhong Zhang       Create
****************************************************************************/

#pragma once
#include <algorithm>
#include <vector>

#include "gtoc9_problem.h"
#include "model.h"

using namespace std;

/****************************************************************************
* Struct       : DP_info_struct
* Description  : struct used in dynamic programming
****************************************************************************/
struct DP_info_struct
{
	double opt_index;
	int last_iter; // iter refers to the number of a specific mission
	double m0;
	double end_epcoh;
	double T_single_mission;
	double dv_single_mission;
};

/****************************************************************************
* Struct       : DP_info_mission
* Description  : struct used in dynamic programming
****************************************************************************/
struct DP_info_mission
{
	double opt_index;
	int iter;
	int end_epoch_iter;  // iter refers to the number of a specific mission
	int current_epoch_iter;  // iter refers to the number of a specific mission
	double epcoh;

	DP_info_mission()
	{
		opt_index = 1.0e10;
		iter = -1;
		current_epoch_iter = -1;
		end_epoch_iter = -1;
		epcoh = -1.0;
	}
};

/****************************************************************************
* Struct       : DP_info_struct_single_mission
* Description  : struct used in dynamic programming
****************************************************************************/
struct DP_info_struct_single_mission
{
	double m0;
	std::vector<double> T_single_mission;
	std::vector<double> dv_single_mission;
};

/****************************************************************************
* Struct       : decision_sequence
* Description  : 
****************************************************************************/
struct decision_sequence
{
	double m;
	int decision[Mission_number];
};

/****************************************************************************
* Function     : DP_optimization_single_mission
* Description  :
*                input: 
*					sequence: debris sequence
*					start_epoch: start epoch (from 0)
*                   end_epoch: end epoch (from 0)
*                ouput: 
*					minimum mass
*					T_single_mission: time for a single task 
****************************************************************************/
vector<DP_info_struct_single_mission> DP_optimization_single_mission(const vector<int>& sequence, double start_epoch, double end_epoch, std::vector<double>& T_single_mission);

/****************************************************************************
* Function     : DP_optimization_single_mission_min_m0
* Description  :
*                input:
*					sequence: debris sequence
*					start_epoch: start epoch (from 0)
*                   end_epoch: end epoch (from 0)
*                ouput:
*					minimum mass
****************************************************************************/
DP_info_struct_single_mission DP_optimization_single_mission_min_m0(const vector<int>& sequence, double start_epoch, double end_epoch);

/****************************************************************************
* Function     : DP_optimization_all_mission
* Description  :
*                input:
*					sequence: debris sequence
*					start_epoch: start epoch (from 0)
*                   end_epoch: end epoch (from 0)
*					opt_min_top: default is 50.0, it is dsigned to accelerate the calculation
*                ouput:
*					T_all: time of all mission
*					dv_all: dv of all mission
****************************************************************************/
double DP_optimization_all_mission(const vector<vector<int>>& sequence, std::vector< std::vector<double>>& T_all, std::vector< std::vector<double>>& dv_all, double
                                   opt_min_top = 50.0);
