/****************************************************************************
* Copyright (C), 2020-2031 清华大学航天航空学院动力学与控制实验室
* 作者: 张众 zhong-zh19@mails.tsinghua.edu.cn
*            539977562@qq.com
* 文件名: main_test_gtoc9.cpp
* 内容简述：使用gtoc9测试蚁群算法文件
*
* 文件历史：
* 版本号     日期         作者       说明
* 01       2021-10-01    张众      创建文件
****************************************************************************/



//数据集debris_data
//0碎片编号 1* 2轨道半长轴(m) 3偏心率 4倾角(rad) 
//5升交点赤经(rad) 6近地点幅角(rad) 7真近点角(rad)
//double debris_data[123][8];   //全局碎片信息

/****************************************************************************
* 函数名   : load_input()
* 功  能   : 读取空间碎片信息,并做预处理
****************************************************************************/

/****************************************************************************
* 函数名   : load_input()
* 功  能   : 根据几次速度增量 计算mp，和初始质量
****************************************************************************/
#include "ACO_gtoc9.h"

vector<vector<int>> divide_sequence(const std::vector<int>& X, int* removal_num_each)
{
	vector<vector<int>> sequence(TreeNum);
	int counter = 0;
	for(int mission =0; mission<TreeNum;mission++)
	{
		for(int i =0; i< removal_num_each[mission];i++)
		{
			sequence[mission].push_back(X[counter]);
			counter++;
		}
	}
	return sequence;
}

//参数：之前的路径;
//      中间过程时间
//      之前的速度增量
double HeurFun_gtoc9(const std::vector<vector<int>>& X_0, const int X_1, void* f_data)
{
	//首先确认没有重复
	bool notin = true;
	for (auto i : X_0)
	{
		for(auto j:i)
		{
			if (X_1 == j)
			{
				notin = false;
				return -1.0;
			}
		}

	}

	auto* para = (Optimization_index*)f_data; //转换为double

	int next_mission = para->expand_mission_ID;
	//if (X_0.size() == 1)
	//{
	//	para[0] = 0.0; para[1] = 0.0; para[500] = 5000.0;
	//}

	//参数中0~499 为ts tf， 当新任务开始时，放入的为tf，tf（第一个点放入0 0）
	//参数中500~999 为dv（每次新发射的为5000）
	//设访问空间碎片的下标编号为n，则ts，tf每个点的时刻为2n,2n+1 （如编号为1，代表第二个点，编号为2,3）
	//每个碎片的dv为500 + n (第一个点，n=0，编号为500)

	double t_now = para->t_mission[next_mission][1];
	double ts, tf;
	double dv = estimate_dv(X_0[next_mission].back(), X_1, t_now, ts, tf);

	if (dv > 800.0) return -1.0;

	double start_epoch, end_epoch;
	start_epoch = para->t_mission[next_mission][0];
	if (next_mission == TreeNum - 1) end_epoch = 2947.0;
	else end_epoch = para->t_mission[next_mission + 1][0] - 35.0;

	double temp_end = (X_0[next_mission].size()) * 23.5 + start_epoch;
	if (temp_end < end_epoch) end_epoch = temp_end;

	vector<int> x_temp = X_0[next_mission];
	x_temp.push_back(X_1);
	auto DP_result = DP_optimization_single_mission_min_m0(x_temp, start_epoch, end_epoch);

	vector<double>& time_sequence = DP_result.T_single_mission;
	vector<double>& dv_sequence = DP_result.dv_single_mission;
	if (time_sequence.size() > 0)
	{
		int i = next_mission;
		{
			//op_index_.t_mission[i][0] = time_sequence[0];
			para->t_mission[i][1] = time_sequence.back();
			para->dv_mission[i] = accumulate(dv_sequence.begin(), dv_sequence.end(), 0.0);
			para->m_mission[i] = DP_result.m0;
			if (DP_result.m0 > 7000.0) return -1.0;

		}
	}

	return 1. / dv_sequence.back();
}



/****************************************************************************
* 函数名   : Obj_gtoc9(const std::vector<int>& X, void* f_data)
* 功  能   : 计算总指标
****************************************************************************/
double Obj_gtoc9_aco(const std::vector<vector<int>>& X, void* f_data)
{

	auto* para = (opt_struct*)f_data;

	vector<vector<double>> T_sequence, dv_missions;
	double total_socre = DP_optimization_total(X, T_sequence, dv_missions);

	int total_num = 0;
	for (int mission = 0; mission < TreeNum; mission++)
	{
		total_num += X[mission].size();
	}

	if(total_num < DebrisNum)
	{
		total_socre += (DebrisNum - total_num) * 550.0;
	}

	para->dv.swap(dv_missions);
	para->t.swap(T_sequence);
	return total_socre;
}



//double localsearch_gtoc9_next(std::vector<vector<int>>& X, void* f_data)
//{
//	auto x0 = X;
//	auto para0 = *(opt_struct*)f_data;
//
//	double dv_max = 2500.0;
//
//	auto* para = (opt_struct*)f_data;
//	double  optimin = para->optimin;
//	auto& dv = para->dv;
//	auto& t = para->t;
//
//	vector<int> visited(DebrisNum,0);
//	for(auto& i: x0)
//	{
//		for(auto j : i)
//		{
//			visited[j] += 1;
//		}
//	}
//
//	//add
//#pragma omp parallel for schedule(dynamic)
//	for (int next_derbis = 0; next_derbis < DebrisNum; next_derbis++)
//	{
//		if(visited[next_derbis] == 0)
//		{
//			for (int i = 0; i < x0.size(); i++)
//			{
//				for (int j = 0; j < x0[i].size()+1; j++)
//				{
//					bool ifswap = true;
//
//					if (j != x0[i].size())
//					{
//						double t_now = para0.t[i][j * 2];
//						double ts, tf;
//						double dv_temp = estimate_dv(next_derbis, x0[i][j], t_now, ts, tf);
//						if (dv_temp > dv_max * 5.0) ifswap = false;
//					}
//
//					if(j != 0)
//					{
//						double t_now = para0.t[i][(j-1) * 2];
//						double ts, tf;
//						double dv_temp = estimate_dv(x0[i][j-1], next_derbis,  t_now, ts, tf);
//						if (dv_temp > dv_max * 5.0) ifswap = false;
//					}
//
//					if(ifswap)
//					{
//						auto x_temp = x0;
//						x_temp[i].insert(x_temp[i].begin() + j, next_derbis);
//						opt_struct a;
//						double opti_temp = Obj_gtoc9_aco(x_temp, &a);
//#pragma omp critical
//                    {
//						if (opti_temp < optimin)
//						{
//							optimin = opti_temp;
//							X.swap(x_temp);
//							dv.swap(a.dv);
//							t.swap(a.t);
//						}
//					}
//
//					}
//				}
//			}
//
//		}
//	}
//
//
//	int counter_num = 0;
//	for(int i =0; i<x0.size();i++)
//	{
//		counter_num += x0[i].size();
//	}
//	if(counter_num != DebrisNum)
//	{
//		return optimin;
//	}
//
//
//
//	//insert
//	//select one debris
//	pair<double, int> pa[123 - TreeNum];
//	int counter = 0;
//	for (int i = 0; i < dv.size(); i++)
//	{
//		for (int j = 0; j < dv[i].size(); j++)
//		{
//			pa[counter] = make_pair(dv[i][j], i * 1000 + j);
//			counter++;
//		}
//	}
//
//	sort(pa, pa + 123 - TreeNum, [](const pair<double, int>& a, const pair<double, int>& b) {return a.first > b.first; });
//
//	auto seque = x0;
//	for (int i = 0; i < x0.size(); i++)
//	{
//		for (int j = 0; j < x0[i].size(); j++)
//		{
//			seque[i][j] = 0;
//		}
//	}
//
//	for (int i = 0; i < 8; i++)
//	{
//		int temp = pa[i].second;
//		seque[temp / 1000][temp % 1000] += 1;
//		seque[temp / 1000][temp % 1000 + 1] += 1;
//	}
//
//	//swap
//#pragma omp parallel for schedule(dynamic)
//	for (int i = 0; i < x0.size(); i++)
//	{
//		for (int j = 0; j < x0[i].size(); j++)
//		{
//			if (seque[i][j] > 0)
//			{
//				for (int k = 0; k < x0.size(); k++)
//				{
//					for (int m = 0; m < x0[k].size(); m++)
//					{
//						bool ifcanswap = true;
//
//						if(j!= x0[i].size()-1)
//						{
//							int next_derbis = x0[k][m];
//							double t_now = para0.t[i][j * 2];
//							double ts, tf;
//							double dv_temp = estimate_dv(next_derbis, x0[i][j+1], t_now, ts, tf);
//							if (dv_temp > dv_max) ifcanswap = false;
//						}
//						if(m != x0[k].size() - 1)
//						{
//							int next_derbis = x0[i][j];
//							double t_now = para0.t[k][m * 2];
//							double ts, tf;
//							double dv_temp = estimate_dv(next_derbis, x0[k][m+1], t_now, ts, tf);
//							if (dv_temp > dv_max) ifcanswap = false;
//
//						}
//						if (m == x0[k].size() - 1 && j == x0[i].size() - 1)
//						{
//							int next_derbis = x0[i][j];
//							double t_now = para0.t[k][(m-1) * 2];
//							double ts, tf;
//							double dv_temp = estimate_dv(x0[k][m - 1], next_derbis , t_now, ts, tf);
//							if (dv_temp > dv_max) ifcanswap = false;
//
//							next_derbis = x0[k][m];
//							t_now = para0.t[i][(j - 1) * 2];
//							ts, tf;
//							dv_temp = estimate_dv(x0[i][j - 1], next_derbis, t_now, ts, tf);
//							if (dv_temp > dv_max) ifcanswap = false;
//						}
//
//						if(ifcanswap && !(i==k && j==m))
//						{
//							auto x_temp = x0;
//							x_temp[i][j] = x0[k][m];
//							x_temp[k][m] = x0[i][j];
//							opt_struct a;
//							double opti_temp = Obj_gtoc9_aco(x_temp, &a);
//#pragma omp critical
//							{
//							if (opti_temp < optimin)
//							{
//								optimin = opti_temp;
//								X.swap(x_temp);
//								dv.swap(a.dv);
//								t.swap(a.t);
//							}
//							}
//						}
//
//
//					}
//
//				}
//			}
//
//
//		}
//	}
//
//	//insert i = k
//#pragma omp parallel for schedule(dynamic)
//	for (int i = 0; i < x0.size(); i++)
//	{
//		for (int j = 0; j < x0[i].size(); j++)
//		{
//			//insert
//			if (seque[i][j] > 0)
//			{
//				int k = i;
//				//for (int k = 0; k < x0.size(); k++)
//				{
//					for (int m = 0; m < x0[k].size(); m++)
//					{
//						bool ifcanswap = true;
//						auto x_temp = x0;
//						int x_insert = x0[i][j];
//						x_temp[i].erase(x_temp[i].begin() + j);
//						x_temp[k].insert(x_temp[k].begin() + m, x_insert);
//
//						if( m!=0)
//						{
//							int next_derbis = x_temp[k][m];
//							double t_now = para0.t[k][((m-1)) * 2];
//							double ts, tf;
//							double dv_temp = estimate_dv(x_temp[k][m - 1], next_derbis,  t_now, ts, tf);
//							if (dv_temp > dv_max) ifcanswap = false;
//
//							//if (m != x_temp[k].size() - 1)
//							//{
//							//	next_derbis = x_temp[k][m];
//							//	t_now = para0.t[k][((m)) * 2];
//							//	ts, tf;
//							//	dv_temp = estimate_dv(next_derbis, x_temp[k][m + 1], t_now, ts, tf);
//							//	if (dv_temp > dv_max) ifcanswap = false;
//							//}
//
//						}
//						else
//						{
//							int next_derbis = x_temp[k][m];
//							double t_now = para0.t[k][(m) * 2];
//							double ts, tf;
//							double dv_temp = estimate_dv(next_derbis, x_temp[k][m + 1], t_now, ts, tf);
//							if (dv_temp > dv_max) ifcanswap = false;
//						}
//
//						if(ifcanswap && (k!=m))
//						{
//							opt_struct a;
//							double opti_temp = Obj_gtoc9_aco(x_temp, &a);
//#pragma omp critical  
//							{
//								if (opti_temp < optimin)
//								{
//									optimin = opti_temp;
//									X.swap(x_temp);
//									dv.swap(a.dv);
//									t.swap(a.t);
//								}
//							}
//						}
//
//
//					}
//
//				}
//			}
//
//		}
//	}
//
//	//insert i != k
//#pragma omp parallel for schedule(dynamic)
//	for (int i = 0; i < x0.size(); i++)
//	{
//		for (int j = 0; j < x0[i].size(); j++)
//		{
//			//insert
//			if (seque[i][j] > 0)
//			{
//				for (int k = 0; k < x0.size(); k++)
//				{
//					if(i == k) continue;
//
//					for (int m = 0; m < x0[k].size()+1; m++)
//					{
//						bool ifcanswap = true;
//						auto x_temp = x0;
//						int x_insert = x0[i][j];
//						x_temp[i].erase(x_temp[i].begin() + j);
//						x_temp[k].insert(x_temp[k].begin() + m, x_insert);
//
//						if (m != 0)
//						{
//							int next_derbis = x_temp[k][m];
//							double t_now = para0.t[k][((m - 1)) * 2];
//							double ts, tf;
//							double dv_temp = estimate_dv(x_temp[k][m - 1], next_derbis, t_now, ts, tf);
//							if (dv_temp > dv_max) ifcanswap = false;
//
//							//if(m != x_temp[k].size()-1)
//							//{
//							//	next_derbis = x_temp[k][m];
//							//	t_now = para0.t[k][((m)) * 2];
//							//	ts, tf;
//							//	dv_temp = estimate_dv(next_derbis, x_temp[k][m + 1], t_now, ts, tf);
//							//	if (dv_temp > dv_max) ifcanswap = false;
//							//}
//
//						}
//						else
//						{
//							int next_derbis = x_temp[k][m];
//							double t_now = para0.t[k][(m) * 2];
//							double ts, tf;
//							double dv_temp = estimate_dv(next_derbis, x_temp[k][m + 1], t_now, ts, tf);
//							if (dv_temp > dv_max) ifcanswap = false;
//						}
//
//						if (ifcanswap && (k != m))
//						{
//							opt_struct a;
//							double opti_temp = Obj_gtoc9_aco(x_temp, &a);
//#pragma omp critical  
//							{
//								if (opti_temp < optimin)
//								{
//									optimin = opti_temp;
//									X.swap(x_temp);
//									dv.swap(a.dv);
//									t.swap(a.t);
//								}
//							}
//						}
//
//
//					}
//
//				}
//			}
//
//		}
//	}
//
//	return optimin;
//}



double localsearch_gtoc9_next_single(std::vector<vector<int>>& X, void* f_data)
{
	auto x0 = X;
	auto para0 = *(opt_struct*)f_data;

	double dv_max = 2500.0;

	auto* para = (opt_struct*)f_data;
	double  optimin = para->optimin;
	auto& dv = para->dv;
	auto& t = para->t;

	vector<int> visited(DebrisNum, 0);
	for (auto& i : x0)
	{
		for (auto j : i)
		{
			visited[j] += 1;
		}
	}

	//add
//#pragma omp parallel for schedule(dynamic)
	for (int next_derbis = 0; next_derbis < DebrisNum; next_derbis++)
	{
		if (visited[next_derbis] == 0)
		{
			for (int i = 0; i < x0.size(); i++)
			{
				for (int j = 0; j < x0[i].size() + 1; j++)
				{
					bool ifswap = true;

					if (j != x0[i].size())
					{
						double t_now = para0.t[i][j * 2];
						double ts, tf;
						double dv_temp = estimate_dv(next_derbis, x0[i][j], t_now, ts, tf);
						if (dv_temp > dv_max * 5.0) ifswap = false;
					}

					if (j != 0)
					{
						double t_now = para0.t[i][(j - 1) * 2];
						double ts, tf;
						double dv_temp = estimate_dv(x0[i][j - 1], next_derbis, t_now, ts, tf);
						if (dv_temp > dv_max * 5.0) ifswap = false;
					}

					if (ifswap)
					{
						auto x_temp = x0;
						x_temp[i].insert(x_temp[i].begin() + j, next_derbis);
						opt_struct a;
						double opti_temp = Obj_gtoc9_aco(x_temp, &a);
//#pragma omp critical
						{
							if (opti_temp < optimin)
							{
								optimin = opti_temp;
								X.swap(x_temp);
								dv.swap(a.dv);
								t.swap(a.t);
							}
						}

					}
				}
			}

		}
	}


	int counter_num = 0;
	for (int i = 0; i < x0.size(); i++)
	{
		counter_num += x0[i].size();
	}
	if (counter_num != DebrisNum)
	{
		return optimin;
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

	for (int i = 0; i < 8; i++)
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
							opt_struct a;
							double opti_temp = Obj_gtoc9_aco(x_temp, &a);
//#pragma omp critical
							{
								if (opti_temp < optimin)
								{
									optimin = opti_temp;
									X.swap(x_temp);
									dv.swap(a.dv);
									t.swap(a.t);
								}
							}
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
				int k = i;
				//for (int k = 0; k < x0.size(); k++)
				{
					for (int m = 0; m < x0[k].size(); m++)
					{
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

							//if (m != x_temp[k].size() - 1)
							//{
							//	next_derbis = x_temp[k][m];
							//	t_now = para0.t[k][((m)) * 2];
							//	ts, tf;
							//	dv_temp = estimate_dv(next_derbis, x_temp[k][m + 1], t_now, ts, tf);
							//	if (dv_temp > dv_max) ifcanswap = false;
							//}

						}
						else
						{
							int next_derbis = x_temp[k][m];
							double t_now = para0.t[k][(m) * 2];
							double ts, tf;
							double dv_temp = estimate_dv(next_derbis, x_temp[k][m + 1], t_now, ts, tf);
							if (dv_temp > dv_max) ifcanswap = false;
						}

						if (ifcanswap && (k != m))
						{
							opt_struct a;
							double opti_temp = Obj_gtoc9_aco(x_temp, &a);
//#pragma omp critical  
							{
								if (opti_temp < optimin)
								{
									optimin = opti_temp;
									X.swap(x_temp);
									dv.swap(a.dv);
									t.swap(a.t);
								}
							}
						}


					}

				}
			}

		}
	}

	//insert i != k
//#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < x0.size(); i++)
	{
		for (int j = 0; j < x0[i].size(); j++)
		{
			//insert
			if (seque[i][j] > 0)
			{
				for (int k = 0; k < x0.size(); k++)
				{
					if (i == k) continue;

					for (int m = 0; m < x0[k].size() + 1; m++)
					{
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

							//if(m != x_temp[k].size()-1)
							//{
							//	next_derbis = x_temp[k][m];
							//	t_now = para0.t[k][((m)) * 2];
							//	ts, tf;
							//	dv_temp = estimate_dv(next_derbis, x_temp[k][m + 1], t_now, ts, tf);
							//	if (dv_temp > dv_max) ifcanswap = false;
							//}

						}
						else
						{
							int next_derbis = x_temp[k][m];
							double t_now = para0.t[k][(m) * 2];
							double ts, tf;
							double dv_temp = estimate_dv(next_derbis, x_temp[k][m + 1], t_now, ts, tf);
							if (dv_temp > dv_max) ifcanswap = false;
						}

						if (ifcanswap && (k != m))
						{
							opt_struct a;
							double opti_temp = Obj_gtoc9_aco(x_temp, &a);
#pragma omp critical  
							{
								if (opti_temp < optimin)
								{
									optimin = opti_temp;
									X.swap(x_temp);
									dv.swap(a.dv);
									t.swap(a.t);
								}
							}
						}


					}

				}
			}

		}
	}

	return optimin;
}
//double localsearch_gtoc9(std::vector<int>& X, void* f_data)
//{
//	std::vector<int> x0 = X;
//
//	double* para = (double*)f_data;
//	double  optimin = para[0];
//
//	double para_temp[1000];
//	double para_opt[1000];
//	//2-opt
//	for (int i = 0; i < x0.size(); i++)
//	{
//		for (int j = i + 1; j < x0.size(); j++)
//		{
//			std::vector<int> x_temp = x0;
//			for (int k = i; k <= j; k++)
//			{
//				x_temp[k] = x0[j + i - k];
//			}
//			double opti_temp = obj_by_sequence(x_temp, para_temp);
//			if (opti_temp < optimin)
//			{
//				optimin = opti_temp;
//				X = x_temp;
//				memcpy(para_temp, para_opt, 1000 * sizeof(double));
//			}
//		}
//	}
//
//	//swap
//	for (int i = 0; i < x0.size(); i++)
//	{
//		for (int j = i + 1; j < x0.size(); j++)
//		{
//			std::vector<int> x_temp = x0;
//			x_temp[i] = x0[j];
//			x_temp[j] = x0[i];
//			double opti_temp = obj_by_sequence(x_temp, para_temp);
//			if (opti_temp < optimin)
//			{
//				optimin = opti_temp;
//				X = x_temp;
//				memcpy(para_temp, para_opt, 1000 * sizeof(double));
//			}
//		}
//	}
//
//	insert
//	for (int i = 0; i < x0.size(); i++)
//	{
//		for (int j = 0; j < x0.size(); j++)
//		{
//			std::vector<int> x_temp = x0;
//			int x_insert = x_temp[i];
//			x_temp.erase(x_temp.begin() + i);
//			x_temp.insert(x_temp.begin() + j, x_insert);
//
//			double opti_temp = obj_by_sequence(x_temp, para_temp);
//			if (opti_temp < optimin)
//			{
//				optimin = opti_temp;
//				X = x_temp;
//				memcpy(para_temp, para_opt, 1000 * sizeof(double));
//			}
//		}
//	}
//
//	f_data = para_opt;
//	return optimin;
//}

double x[123] = { 23,55,79,113,25,20,27,117,121,50,95,102,38,97,
 19,115,41,26,45,82,47,85,7,2,11,77,
	72,107,61,10,28,3,64,66,31,90,73,87,57,35,69,65,8,43,71,4,29,
	108,24,104,119,22,75,63,112,37,32,114,
	84,59,98,1,40,51,36,67,62,99,54,122,76,15,
	101,48,53,5,12,39,58,13,60,74,
	49,9,70,93,105,46,88,118,18,91,
	86,34,100,30,92,6,110,96,81,
	33,68,116,106,14,52,120,80,16,94,83,89,
	44,111,56,78,0,17,109,103,42,21
};

//void main()
//{
//
//	//std::vector<int> x_vector(x,x+123);
//	//double para[1000];
//	//std::cout<< obj_by_sequence(x_vector, para);
//
//	ACO aco(Obj_gtoc9, HeurFun_gtoc9, localsearch_gtoc9, 122);
//
//	aco.Run();
//}