#pragma once
#include <vector>

#include "Constant.h"

#include "DP.h"

using namespace std;

double mp_calc(const std::vector<double>& dv, double& m0);
std::vector<double> divide_missions(const vector<vector<int>>& debris_sequence);
void map_input_gtoc9(vector<double>& X);
void inverse_map_input_gtoc9(vector<double>& X);
double penalty_gtoc9(const vector<vector<double>>& T);
/****************************************************************************
* Struct       : Opt_info_gtoc9
* Discription  : Used to save sequence info
****************************************************************************/
struct Opt_info_gtoc9
{
	int mission;
	double end_peoch;
	int last_debris_num;
	vector<vector<double>> time_sequence;
	vector<vector<int>> debris_squence;
};


struct opt_struct
{
	double  optimin;
	double opt_each_mission[TreeNum];
	std::vector<vector<double>> t;
	std::vector<vector<double>> dv;
};

double Obj_gtoc9(const std::vector<double>& X, std::vector<double>& grad, void* f_data);
double Obj_gtoc9_single_mission_single_debris(const std::vector<double>& X, std::vector<double>& grad, void* f_data);
double Obj_gtoc9_single_mission(const std::vector<double>& X, std::vector<double>& grad, void* f_data);
vector<double> out_gtoc9_index_m0(const std::vector<double>& X, void* f_data);
vector<double> estimate_t_sequence(const vector<vector<int>>& debris_sequence);
//vector<double> estimate_t_sequence_random(const vector<vector<int>>& debris_sequence);
//void optimize_single_mission(vector<vector<double>>& T_sequence, Opt_info_gtoc9& opt_info, int N);
//void optimize_single_mission_ts_tf(vector<vector<double>>& T_sequence, Opt_info_gtoc9& opt_info, int N);
//
//void single_mission(const vector<double>& mission_start_epoch, const Opt_info_gtoc9& opt_info, vector< vector<double>>& T_sequence);
//double time_optimization_NLOPT(vector<double>& T_opt, const vector<vector<int>>& debris_sequence);