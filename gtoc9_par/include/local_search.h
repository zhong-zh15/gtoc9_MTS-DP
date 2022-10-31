/****************************************************************************
* Copyright (C), 2020-2031 Tsinghua University, School of Aerospace Engineering, LAD
* Author: Zhong Zhang
				zhong-zh19@mails.tsinghua.edu.cn
*               539977562@qq.com
* File: local_search.h
* Description: call the optimization function to find the optimal solution
*
* Log:
*Version      Date        Author           Description
* 01        2022-03-15    Zhong Zhang       Create
****************************************************************************/

#ifndef LOCALSEARCH
#define LOCALSEARCH
#include<iostream>
#include<string>
#include <unordered_map>
#include<vector>
#include <unordered_set>

#include "gtoc9_problem.h"

using namespace std;

/****************************************************************************
* Struct       : each_mis_neighborhood
* Description  : mission_neighborhood is a struct that represents a neighborhood of a given mission
****************************************************************************/
struct each_mis_neighborhood
{
	vector<pair<string, int>> neighborhood_each;
};

/****************************************************************************
* Function     : string2sequence
* Description  : return all mission debris sequences of a given string	
*                input:
*					x: a given string
*                ouput:
*					return: all mission debris sequences of a given string	
****************************************************************************/
std::vector<std::vector<int>> string2sequence(const string& x);

/****************************************************************************
* Function     : string2sequence
* Description  : return single mission debris sequence of a given string and the mission id
*                input:
*					x: a given string
*					mission_id:  mission id
*                ouput:
*					return: single mission debris sequence of a given string and the mission id
****************************************************************************/
std::vector<int> string2sequence(const string& x, int& mission_id);

/****************************************************************************
* Function     : sequence2string
* Description  : tranfer a given all mission sequence to a string
*                input:
*					x: a given all mission sequence
*                ouput:
*					return: a string
****************************************************************************/
string sequence2string(const std::vector<std::vector<int>>& x);

/****************************************************************************
* Function     : sequence2string
* Description  : tranfer a given single mission sequence and the mission ID to a string
*                input:
*					x: a given single mission sequence
*					mission: mission ID
*                ouput:
*					return: a string
****************************************************************************/
string sequence2string(const std::vector<int>& x, int mission);

/****************************************************************************
* Function     : local_search_1layer
* Description  : local search algorithm for gtoc9 problem called in multitree_beam.cpp-Localsearch function
*                compute the best solution for each mission sequence in its neighborhood
*                input:
*					neighborhood: a database for all computed mission neighborhoods
*					X_all: all mission sequences
*                   a_all: the corrsponding info (time sequence and delta v sequence) to a_all
*					opt_min: the maximum cost in X_all (a parameter used in DP to accelerate calculation)
*                ouput:
*				    X_all: the best solution for each mission sequence in its neighborhood
*                   a_all: the corrsponding info (time sequence and delta v sequence) to a_all
*					X_all_end: the local optimal sequence (if one sequence cannot find a better solution in its neighborhood, move it to X_all_end)
*                   a_all_end: the corrsponding info (time sequence and delta v sequence) to X_all_end
****************************************************************************/
void local_search_1layer(unordered_map<string, double>& neighborhood, vector<vector<vector<int>>>& X_all,
	vector<opt_struct>& a_all, vector<vector<vector<int>>>& X_all_end,
	vector<opt_struct>& a_all_end, double opt_min);

/****************************************************************************
* Function     : localsearch_gtoc9_MTS_pool
* Description  : compute the neighborhood for single mission sequence
*                input:
*					X: mission sequence
*                   f_data: the corrsponding info (time sequence and delta v sequence) to X
*					neighbor:  a database for all computed mission neighborhoods
*                ouput:
*				    return value: the neighborhood for each single mission sequence (in string)
****************************************************************************/
vector<string> localsearch_gtoc9_MTS_pool(std::vector<vector<int>>& X, opt_struct& f_data, vector<each_mis_neighborhood>& neighbor);

/****************************************************************************
* Function     : localsearch_gtoc9_MTS_changetime_pool
* Description  :
*                compute the neighborhood for all mission sequence
*                input:
*					X: mission sequence
*                   f_data: the corrsponding info (time sequence and delta v sequence) to X
*                ouput:
*				    return value: the neighborhood for all mission sequence (in string)
****************************************************************************/
vector<string> localsearch_gtoc9_MTS_changetime_pool(std::vector<vector<int>>& X, opt_struct& f_data);

#endif