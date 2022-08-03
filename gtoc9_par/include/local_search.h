#ifndef LOCALSEARCH
#define LOCALSEARCH
#include<iostream>
#include<string>
#include <unordered_map>
#include<vector>
#include <unordered_set>

//#include "ACO_gtoc9.h"

#include "gtoc9_problem.h"

using namespace std;

struct each_mis_neighborhood
{
	vector<pair<string, int>> neighborhood_each;
};

void local_search_1layer(unordered_map<string, double>& neighborhood, vector<vector<vector<int>>>& X_all,
	vector<opt_struct>& a_all, vector<vector<vector<int>>>& X_all_end,
	vector<opt_struct>& a_all_end, double opt_min);

vector<string> localsearch_gtoc9_MTS_pool(std::vector<vector<int>>& X, opt_struct& f_data, vector<each_mis_neighborhood>& neighbor);

vector<string> localsearch_gtoc9_MTS_changetime_pool(std::vector<vector<int>>& X, opt_struct& f_data);
std::vector<std::vector<int>> string2sequence(const string &x);

std::vector<int> string2sequence(const string& x, int &mission_id);

string sequence2string(const std::vector<std::vector<int>>& x);
string sequence2string(const std::vector<int>& x, int mission);
#endif
