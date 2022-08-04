#pragma once
#include <vector>

#include "Constant.h"

#include "DP.h"

using namespace std;

double mp_calc(const std::vector<double>& dv, double& m0);

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

/****************************************************************************
* Struct       : opt_struct
* Discription  : save each mission's info
****************************************************************************/
struct opt_struct
{
	double  optimin;
	double opt_each_mission[TreeNum];
	std::vector<vector<double>> t;
	std::vector<vector<double>> dv;
};