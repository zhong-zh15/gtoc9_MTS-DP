#pragma once

#include <fstream>
#include <numeric>

#include "OrbitFun.h"
#include "Constant.h"
#include "model.h"
#include "aco.h"
#include "gtoc9_problem.h"
#include "problem_struct.h"

//double localsearch_gtoc9_next(std::vector<vector<int>>& X, void* f_data);
double localsearch_gtoc9_next_single(std::vector<vector<int>>& X, void* f_data);
double HeurFun_gtoc9(const std::vector<vector<int>>& X_0, const int X_1, void* f_data);
struct opt_struct
{
	double  optimin;
	double opt_each_mission[TreeNum];
	std::vector<vector<double>> t;
	std::vector<vector<double>> dv;
};

double Obj_gtoc9_aco(const std::vector<vector<int>>& X, void* f_data);
//double localsearch_gtoc9_next(std::vector<vector<int>>& X, void* f_data);