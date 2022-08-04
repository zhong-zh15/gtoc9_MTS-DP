#pragma once
#include <algorithm>
#include <vector>

#include "gtoc9_problem.h"
#include "model.h"

using namespace std;

struct DP_info_struct
{
	double opt_index;
	int last_iter;
	double m0;
	double end_epcoh;
	double T_single_mission;
	double dv_single_mission;
};


struct DP_info_mission
{
	double opt_index;
	int iter;
	int end_epoch_iter;
	int current_epoch_iter;
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

struct DP_info_struct_single_mission
{
	double m0;
	std::vector<double> T_single_mission;
	std::vector<double> dv_single_mission;
};

struct decision_sequence
{
	double m;
	int decision[Mission_number];
};


/****************************************************************************
* Function     : DP_optimization
* Discription  :
*                input: sequence()
*				        start_epoch
*                       end_epoch ()
*                ouput: T_single_mission
****************************************************************************/
vector<DP_info_struct_single_mission> DP_optimization_single_mission(const vector<int>& sequence, double start_epoch, double end_epoch, std::vector<double>& T_single_mission);

DP_info_struct_single_mission DP_optimization_single_mission_min_m0(const vector<int>& sequence, double start_epoch, double end_epoch);

double DP_optimization_all_mission(const vector<vector<int>>& sequence, std::vector< std::vector<double>>& T_all, std::vector< std::vector<double>>& dv_all, double
                                   opt_min_top = 50.0);
