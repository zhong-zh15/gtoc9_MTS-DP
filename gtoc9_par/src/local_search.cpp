/****************************************************************************
* Copyright (C), 2020-2031 Tsinghua University, School of Aerospace Engineering, LAD
* Author: Zhong Zhang
				zhong-zh19@mails.tsinghua.edu.cn
*               539977562@qq.com
* File: local_search.cpp
* Description: call the optimization function to find the optimal solution
*
* Log:
*Version      Date        Author           Description
* 01        2022-03-15    Zhong Zhang       Create
****************************************************************************/

#include "local_search.h"

#include <numeric>

/****************************************************************************
* Function     : splitStr
* Description  : split string to int
*                input:
*					current_string: a given string
*					cut: cutted string
*                ouput:
*					arr: int vector
****************************************************************************/
void splitStr(const string& str, vector<int>& arr, const string& cut)
{
	// str: current_string
	// arr: output
	// cut; cut_string
	string::size_type pos1, pos2;
	pos2 = str.find(cut);
	pos1 = 0;
	while (string::npos != pos2)
	{
		arr.push_back(stoi(str.substr(pos1, pos2 - pos1)));
		pos1 = pos2 + cut.size();
		pos2 = str.find(cut, pos1);
	}
	if (pos1 != str.length())
		arr.push_back(stoi(str.substr(pos1)));
}

/****************************************************************************
* Function     : string2sequence
* Description  : return all mission debris sequences of a given string
*                input:
*					x: a given string
*                ouput:
*					return: all mission debris sequences of a given string
****************************************************************************/
std::vector<std::vector<int>> string2sequence(const string& x)
{
	vector<vector<int>> sequence;
	vector<int> all_number;
	string str1 = ",";
	splitStr(x, all_number, str1);
	vector<int> single_sequence;
	int counter = 0;
	while (true)
	{
		int temp = all_number[counter];
		single_sequence.resize(temp);
		for (int i = 0; i < temp; i++)
		{
			counter++;
			single_sequence[i] = all_number[counter];
		}
		counter++;
		sequence.push_back(single_sequence);
		single_sequence.clear();
		if (counter == all_number.size()) break;
	}
	return sequence;
}

/****************************************************************************
* Function     : string2sequence
* Description  : return single mission debris sequence of a given string and the mission id
*                input:
*					x: a given string
*					mission_id: mission id
*                ouput:
*					return: single mission debris sequence of a given string and the mission id
****************************************************************************/
std::vector<int> string2sequence(const string& x, int& mission_id)
{
	vector<int> sequence;
	vector<int> all_number;
	string str1 = ",";
	splitStr(x, all_number, str1);
	int counter = 0;


	mission_id = all_number[counter];
	counter++;

	int temp = all_number[counter];
	sequence.resize(temp);

	for (int i = 0; i < temp; i++)
	{
		counter++;
		sequence[i] = all_number[counter];
	}
	counter++;
	return sequence;
}

/****************************************************************************
* Function     : sequence2string
* Description  : tranfer a given all mission sequence to a string
*                input:
*					x: a given all mission sequence
*                ouput:
*					return: a string
****************************************************************************/
string sequence2string(const std::vector<std::vector<int>>& x)
{
	string string_sequence;

	for (int mission = 0; mission < x.size(); mission++)
	{
		string_sequence += to_string(x[mission].size());
		string_sequence += ",";
		for (int i = 0; i < x[mission].size(); i++)
		{
			string_sequence += to_string(x[mission][i]);
			if (!(mission == x.size() - 1 && i == x[mission].size() - 1)) string_sequence += ",";
		}
	}
	return string_sequence;
}

/****************************************************************************
* Function     : sequence2string
* Description  : tranfer a given single mission sequence and the mission ID to a string
*                input:
*					x: a given single mission sequence
*					mission: mission ID
*                ouput:
*					return: a string
****************************************************************************/
string sequence2string(const std::vector<int>& x, int mission)
{
	string string_sequence;

	string_sequence += to_string(mission);
	string_sequence += ",";

	string_sequence += to_string(x.size());
	string_sequence += ",";
	for (int i = 0; i < x.size(); i++)
	{
		string_sequence += to_string(x[i]);
		if (!(mission == x.size() - 1 && i == x.size() - 1)) string_sequence += ",";
	}
	return string_sequence;
}

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
void local_search_1layer(unordered_map<string, double> &neighborhood, vector<vector<vector<int>>> &X_all, vector<opt_struct> &a_all, vector<vector<vector<int>>> &X_all_end, vector<opt_struct> &a_all_end, double opt_min)
{
	vector<vector<vector<int>>> X_all_next;
	vector<opt_struct> a_all_next;
	
	//compute neighborhood and combine what not in database to one
	unordered_set<string> each_neighborhood_layer;
	vector<vector<string>> all_neighborhood_layer(X_all.size());

	vector<vector<each_mis_neighborhood>> all_neighborhood_layer_singlemission(X_all.size());

	int remove_num = 0;
	auto& sequence_temp = X_all[0];
	for (int i = 0; i < sequence_temp.size(); i++) remove_num += sequence_temp[i].size();

#pragma omp parallel for schedule(dynamic)
	for (int num = 0; num < X_all.size(); num++)
	{
		if (remove_num >= DebrisNum - Last_mission)
		    all_neighborhood_layer[num] = localsearch_gtoc9_MTS_changetime_pool(X_all[num], a_all[num]);
		else
			all_neighborhood_layer[num] = localsearch_gtoc9_MTS_pool(X_all[num], a_all[num], all_neighborhood_layer_singlemission[num]);


	}


	for (int num = 0; num < X_all.size(); num++)
	{
		auto& each_neighborhood = all_neighborhood_layer[num];
		for (int i = 0; i < each_neighborhood.size(); i++)
		{
			if (neighborhood.count(each_neighborhood[i]) == 0)
			{ //not in database, need to compute later
				each_neighborhood_layer.emplace(each_neighborhood[i]);
			}
		}
	}

	//compute cost and insert in database
	vector<string> each_neighborhood_layer_vector; each_neighborhood_layer_vector.reserve(each_neighborhood_layer.size());
	for (auto &iter:each_neighborhood_layer) each_neighborhood_layer_vector.push_back(iter);

	vector<string> string_temp(each_neighborhood_layer_vector.size());
	vector<double> cost_temp(each_neighborhood_layer_vector.size());
	
#pragma omp parallel for schedule(dynamic)
	for (int num = 0; num < each_neighborhood_layer_vector.size(); num++)
	{
		string str1 = each_neighborhood_layer_vector[num];
		double cost_score;
		if (remove_num >= DebrisNum - Last_mission)
		{
			auto sequence = string2sequence(str1);
			vector<vector<double>> t, dv;
			cost_score = DP_optimization_all_mission(sequence, t, dv, opt_min);
		}
		else
		{
			int mission_id = -1;
			auto sequence = string2sequence(str1, mission_id);
			double t_start, t_end;
			t_start =  (2947.0 + 35.0) / TreeNum * mission_id;
			t_end   =  (2947.0 + 35.0) / TreeNum * (mission_id+1.0) -35.0;
			if (mission_id == TreeNum - 1) t_end = 2947.0;
			auto result = DP_optimization_single_mission_min_m0(sequence, t_start, t_end);
			cost_score = 55.0 + 2.0e-6 * (result.m0 - 2000.0) * (result.m0 - 2000.0);
		}
		string_temp[num] = str1;
		cost_temp[num] = cost_score;
	}

	for (int num = 0; num < each_neighborhood_layer_vector.size(); num++)
	{
		neighborhood.emplace(string_temp[num], cost_temp[num]);
	}
	
	//select best one in neighborhood
	if (remove_num >= DebrisNum - Last_mission)
	{
#pragma omp parallel for schedule(dynamic)
		for (int num = 0; num < X_all.size(); num++)
		{
			auto& each_neighborhood = all_neighborhood_layer[num];
			vector<vector<int>> best_sequence;
			double min_cost = 1.0e40;
			for (int i = 0; i < each_neighborhood.size(); i++)
			{
				double cost_temp = neighborhood.find(each_neighborhood[i])->second;
				if (cost_temp < min_cost)
				{
					min_cost = cost_temp;
					best_sequence = string2sequence(each_neighborhood[i]);
				}
			}

			if (min_cost < a_all[num].optimin - 1.0e-6)
			{
				opt_struct a_temp;
				a_temp.optimin = DP_optimization_all_mission(best_sequence, a_temp.t, a_temp.dv, opt_min);
#pragma omp critical
				{
					X_all_next.push_back(best_sequence);
					a_all_next.push_back(a_temp);
				}
			}
			else
			{ //ready to output
				//opt_struct a_temp;
				//a_temp.optimin = DP_optimization_all_mission(best_sequence, a_temp.t, a_temp.dv);
#pragma omp critical
				{
					if (a_all[num].t.size() > 0)
					{
						X_all_end.push_back(X_all[num]);
						a_all_end.push_back(a_all[num]);
					}
				}
			}
		}
	}
	else
	{
#pragma omp parallel for schedule(dynamic)
		for (int num = 0; num < X_all.size(); num++)
		{
			auto& each_neighborhood = all_neighborhood_layer_singlemission[num];
			vector<vector<int>> best_sequence(X_all[num]);
			vector < double > best_cost_eachmission(TreeNum);
			int best_nerghber_iter;
			double min_cost = 1.0e20;
			for (int i = 0; i < each_neighborhood.size(); i++)
			{
				vector<vector<int>> sequence_temp(X_all[num]);
				vector<double> cost_each_mission(TreeNum);
				for (int mission = 0; mission < TreeNum; mission++) cost_each_mission[mission] = a_all[num].opt_each_mission[mission];

				for(int j =0; j< each_neighborhood[i].neighborhood_each.size(); j++)
				{
					auto & single_mission_ = each_neighborhood[i].neighborhood_each[j];
					cost_each_mission[single_mission_.second] = neighborhood.find(single_mission_.first)->second;
					int temp;
					sequence_temp[single_mission_.second] = string2sequence(single_mission_.first, temp);
				}
				double cost_temp = accumulate(cost_each_mission.begin(), cost_each_mission.end(),0.0);
				if (cost_temp < min_cost)
				{
					min_cost = cost_temp;
					best_sequence = sequence_temp;
					best_cost_eachmission = cost_each_mission;
					best_nerghber_iter = i;
				}
			}

			if (min_cost < a_all[num].optimin)
			{
				opt_struct a_temp = a_all[num];
				a_temp.optimin = min_cost;
				for (int i = 0; i < TreeNum; i++) a_temp.opt_each_mission[i] = best_cost_eachmission[i];
				auto& neighbor_mission = each_neighborhood[best_nerghber_iter].neighborhood_each;
				for(int i =0; i< neighbor_mission.size(); i++)
				{
					int mission_current = neighbor_mission[i].second; 
					double t_start, t_end;
					t_start = (2947.0 + 35.0) / TreeNum * mission_current;
					t_end = (2947.0 + 35.0) / TreeNum * (mission_current + 1.0) - 35.0;
					if (mission_current == TreeNum - 1) t_end = 2947.0;
					auto result = DP_optimization_single_mission_min_m0(best_sequence[mission_current], t_start, t_end);
					a_temp.t[mission_current] = result.T_single_mission;
					a_temp.dv[mission_current] = result.dv_single_mission;
				}
#pragma omp critical
				{
					X_all_next.push_back(best_sequence);
					a_all_next.push_back(a_temp);
				}
			}
			else
			{ //ready to output
				//opt_struct a_temp;
				//a_temp.optimin = DP_optimization_all_mission(best_sequence, a_temp.t, a_temp.dv);
#pragma omp critical
				{
					if (a_all[num].t.size() > 0)
					{
						X_all_end.push_back(X_all[num]);
						a_all_end.push_back(a_all[num]);
					}
				}
			}
		}

	}

	X_all.swap(X_all_next);
	a_all.swap(a_all_next);
	X_all_next.clear();
	a_all_next.clear();
}

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
vector<string> localsearch_gtoc9_MTS_pool(std::vector<vector<int>>& X, opt_struct& f_data, vector<each_mis_neighborhood> &neighbor)
{

	
	vector<string> neighborhood_single_mission;
	auto x0 = X;
	auto para0 = f_data;

	double dv_max = 3000.0;

	vector<int> visited(DebrisNum, 0);
	for (auto& i : x0)
	{
		for (auto j : i)
		{
			visited[j] += 1;
		}
	}

	//swap
//#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < x0.size(); i++)
	{
		for (int j = 0; j < x0[i].size(); j++)
		{

				for (int k = 0; k < x0.size(); k++)
				{
					for (int m = 0; m < x0[k].size(); m++)
					{

						if (x0[i].size() < 3 || x0[k].size() < 3) continue;

						bool ifcanswap = true;

						if (j != x0[i].size() - 1)
						{
							int next_derbis = x0[k][m];
							double t_now = para0.t[i][j * 2];
							double ts, tf;
							double dv_temp = estimate_dv(next_derbis, x0[i][j + 1], t_now, ts, tf);
							if (dv_temp > dv_max) ifcanswap = false;
						}
						if (m != x0[k].size() - 1)
						{
							int next_derbis = x0[i][j];
							double t_now = para0.t[k][m * 2];
							double ts, tf;
							double dv_temp = estimate_dv(next_derbis, x0[k][m + 1], t_now, ts, tf);
							if (dv_temp > dv_max) ifcanswap = false;

						}
						if (m == x0[k].size() - 1 && j == x0[i].size() - 1)
						{
							int next_derbis = x0[i][j];
							double t_now = para0.t[k][(m - 1) * 2];
							double ts, tf;
							double dv_temp = estimate_dv(x0[k][m - 1], next_derbis, t_now, ts, tf);
							if (dv_temp > dv_max) ifcanswap = false;

							next_derbis = x0[k][m];
							t_now = para0.t[i][(j - 1) * 2];
							ts, tf;
							dv_temp = estimate_dv(x0[i][j - 1], next_derbis, t_now, ts, tf);
							if (dv_temp > dv_max) ifcanswap = false;
						}

						if (ifcanswap && !(i == k && j == m))
						{
							auto x_temp = x0;
							x_temp[i][j] = x0[k][m];
							x_temp[k][m] = x0[i][j];

							
							{
								each_mis_neighborhood neighborhood_each;
								neighborhood_single_mission.emplace_back(sequence2string(x_temp[i], i));
								neighborhood_each.neighborhood_each.push_back(make_pair(sequence2string(x_temp[i], i), i));
								if (i != k)
								{
									neighborhood_single_mission.emplace_back(sequence2string(x_temp[k], k));
									neighborhood_each.neighborhood_each.push_back(make_pair(sequence2string(x_temp[k], k), k));
								}
								neighbor.push_back(neighborhood_each);
							}
						
						}


					}

				}

		}
	}

	//insert i = k
//#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < x0.size(); i++)
	{
		for (int j = 0; j < x0[i].size(); j++)
		{
			
			
				
				for (int k = 0; k < x0.size(); k++)
				{
					for (int m = 0; m < x0[k].size() + 1; m++)
					{

						if (x0[i].size() < 3 || x0[k].size() < 3) continue;

						if (i == k && m == x0[k].size()) continue;


						bool ifcanswap = true;
						auto x_temp = x0;
						int x_insert = x0[i][j];
						x_temp[i].erase(x_temp[i].begin() + j);
						x_temp[k].insert(x_temp[k].begin() + m, x_insert);

						if (m != 0)
						{
							int next_derbis = x_temp[k][m];
							double t_now = para0.t[k][((m - 1)) * 2];
							double ts, tf;
							double dv_temp = estimate_dv(x_temp[k][m - 1], next_derbis, t_now, ts, tf);
							if (dv_temp > dv_max) ifcanswap = false;



						}
						else
						{
							int next_derbis = x_temp[k][m];
							double t_now = para0.t[k][(m) * 2];
							double ts, tf;
							double dv_temp = estimate_dv(next_derbis, x_temp[k][m + 1], t_now, ts, tf);
							if (dv_temp > dv_max) ifcanswap = false;
						}

						if (ifcanswap && !(i == k && j == m))
						{
							
							{

								each_mis_neighborhood neighborhood_each;
								neighborhood_single_mission.emplace_back(sequence2string(x_temp[i], i));
								neighborhood_each.neighborhood_each.push_back(make_pair(sequence2string(x_temp[i], i), i));
								if (i != k)
								{
									neighborhood_single_mission.emplace_back(sequence2string(x_temp[k], k));
									neighborhood_each.neighborhood_each.push_back(make_pair(sequence2string(x_temp[k], k), k));
								}
								neighbor.push_back(neighborhood_each);
							}
						}


					}

				}
			

		}
	}


	return neighborhood_single_mission;
}

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
vector<string> localsearch_gtoc9_MTS_changetime_pool(std::vector<vector<int>>& X, opt_struct& f_data)
{
	vector<string> local_field;
	
	auto x0 = X;
	auto para0 = f_data;

	if(f_data.t.size() ==0 )
	{
		return local_field;
	}
	else
	{
		if (f_data.t[0].size() ==0)
		{
			return local_field;
		}
	}

	//double dv_max = 1500.0;
	//int range = 30;

	double dv_max = 3000.0;
	int range = 110;

	//auto* para = &f_data;
	double  optimin = f_data.optimin;
	auto& dv = f_data.dv;
	auto& t = f_data.t;
	vector<int> visited(DebrisNum, 0);
	for (auto& i : x0)
	{
		for (auto j : i)
		{
			visited[j] += 1;
		}
	}

	//insert
	//select one debris
	pair<double, int> pa[123 - TreeNum];
	int counter = 0;
	for (int i = 0; i < dv.size(); i++)
	{
		for (int j = 0; j < dv[i].size(); j++)
		{
			pa[counter] = make_pair(dv[i][j], i * 1000 + j);
			counter++;
		}
	}
	sort(pa, pa + 123 - TreeNum, [](const pair<double, int>& a, const pair<double, int>& b) {return a.first > b.first; });

	auto seque = x0;
	for (int i = 0; i < x0.size(); i++)
	{
		for (int j = 0; j < x0[i].size(); j++)
		{
			seque[i][j] = 0;
		}
	}

	for (int i = 0; i < range; i++)
	{
		int temp = pa[i].second;
		seque[temp / 1000][temp % 1000] += 1;
		seque[temp / 1000][temp % 1000 + 1] += 1;
	}


	//swap
//#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < x0.size(); i++)
	{
		for (int j = 0; j < x0[i].size(); j++)
		{
			if (seque[i][j] > 0)
			{
				for (int k = 0; k < x0.size(); k++)
				{
					for (int m = 0; m < x0[k].size(); m++)
					{

						if (x0[i].size() < 3 || x0[k].size() < 3) continue;

						bool ifcanswap = true;

						if (j != x0[i].size() - 1)
						{
							int next_derbis = x0[k][m];
							double t_now = para0.t[i][j * 2];
							double ts, tf;
							double dv_temp = estimate_dv(next_derbis, x0[i][j + 1], t_now, ts, tf);
							if (dv_temp > dv_max) ifcanswap = false;
						}
						if (m != x0[k].size() - 1)
						{
							int next_derbis = x0[i][j];
							double t_now = para0.t[k][m * 2];
							double ts, tf;
							double dv_temp = estimate_dv(next_derbis, x0[k][m + 1], t_now, ts, tf);
							if (dv_temp > dv_max) ifcanswap = false;

						}
						if (m == x0[k].size() - 1 && j == x0[i].size() - 1)
						{
							int next_derbis = x0[i][j];
							double t_now = para0.t[k][(m - 1) * 2];
							double ts, tf;
							double dv_temp = estimate_dv(x0[k][m - 1], next_derbis, t_now, ts, tf);
							if (dv_temp > dv_max) ifcanswap = false;

							next_derbis = x0[k][m];
							t_now = para0.t[i][(j - 1) * 2];
							ts, tf;
							dv_temp = estimate_dv(x0[i][j - 1], next_derbis, t_now, ts, tf);
							if (dv_temp > dv_max) ifcanswap = false;
						}

						if (ifcanswap && !(i == k && j == m))
						{
							auto x_temp = x0;
							x_temp[i][j] = x0[k][m];
							x_temp[k][m] = x0[i][j];
							local_field.emplace_back(sequence2string(x_temp));
						}


					}

				}
			}


		}
	}

	//insert i = k
//#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < x0.size(); i++)
	{
		for (int j = 0; j < x0[i].size(); j++)
		{
			//insert
			if (seque[i][j] > 0)
			{
				//int k = i;
				for (int k = 0; k < x0.size(); k++)
				{
					for (int m = 0; m < x0[k].size() + 1; m++)
					{

						if (x0[i].size() < 3 || x0[k].size() < 3) continue;

						if (i == k && m == x0[k].size()) continue;


						bool ifcanswap = true;
						auto x_temp = x0;
						int x_insert = x0[i][j];
						x_temp[i].erase(x_temp[i].begin() + j);
						x_temp[k].insert(x_temp[k].begin() + m, x_insert);

						if (m != 0)
						{
							int next_derbis = x_temp[k][m];
							double t_now = para0.t[k][((m - 1)) * 2];
							double ts, tf;
							double dv_temp = estimate_dv(x_temp[k][m - 1], next_derbis, t_now, ts, tf);
							if (dv_temp > dv_max) ifcanswap = false;
						}
						else
						{
							int next_derbis = x_temp[k][m];
							double t_now = para0.t[k][(m) * 2];
							double ts, tf;
							double dv_temp = estimate_dv(next_derbis, x_temp[k][m + 1], t_now, ts, tf);
							if (dv_temp > dv_max) ifcanswap = false;
						}

						if (ifcanswap && !(i == k && j == m))
						{
							local_field.emplace_back(sequence2string(x_temp));
						}
					}

				}
			}

		}
	}

	return local_field;
}
