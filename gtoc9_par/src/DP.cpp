/****************************************************************************
* Copyright (C), 2020-2031 Tsinghua University, School of Aerospace Engineering, LAD
* Author: Zhong Zhang
				zhong-zh19@mails.tsinghua.edu.cn
*               539977562@qq.com
* File: DP.cpp
* Description: the specific implementation of dynamic programming on the GTOC9 problem
*
* Log:
*Version      Date        Author           Description
* 01        2022-03-15    Zhong Zhang       Create
****************************************************************************/

#include "DP.h"

#include <set>
#include <cstring>

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
DP_info_struct_single_mission DP_optimization_single_mission_min_m0(const vector<int>& sequence, double start_epoch, double end_epoch)// , std::vector<double>& T_single_mission)
{
	double step_init = 8.0; //days
	int num_step = 10;
	double t_step_tf = 8.0;
	int tf_num = 0;           //discrete num
	double loss_factor = 0.3;
	for (double tf_transfer_temp = 6.0; tf_transfer_temp < 30.0 + 1.0e-1; tf_transfer_temp += t_step_tf) tf_num++;
	int N = sequence.size() - 1;  //decision num

	vector<vector<DP_info_struct>> dp_info(N); //optimum_state
	for (int i = 0; i < N; i++) dp_info[i].reserve(50);

	vector<double> start_T_sequence;
	for (int i = 0; i < num_step; i++)
	{
		start_T_sequence.push_back(start_epoch + i * step_init);
	}

	if (sequence.size() == 1)
	{
		DP_info_struct_single_mission temp;
		temp.m0 = 2000;
		temp.T_single_mission.push_back(start_epoch);
		//temp.dv_single_mission.push_back(0.0);
		return temp;
	}

	//each end epoch
	int now_debris = sequence[0];
	int target_debris = sequence[1];
	int counter = 0;

	//determine x0
	set<double> init_time;
	for (int i = 0; i < start_T_sequence.size(); i++)
	{
		for (int j = 0; j < tf_num; j++)
		{
			init_time.insert(start_T_sequence[i] + 6.0 + j * t_step_tf);
		}
	}
	for (auto iter = init_time.begin(); iter != init_time.end(); ++iter)
	{
		double t_f = *iter;
		double opt_index_min = 1.0e6;
		double T_now, Ts, Tf;
		double dv_min = 1.0e6;
		for (int i = 0; i < start_T_sequence.size(); i++)
		{
			double t_transfer = t_f - start_T_sequence[i] - 6.0;
			if (!(0.0 - 1.0e-6 < t_transfer && t_transfer < 24.0 + 1.0e-6)) continue;
			double t_temp = t_transfer / t_step_tf;
			t_temp = t_temp - round(t_temp);
			if (fabs(t_temp) > 1.0e-3) continue;
			double dv = Dv_ij(now_debris, target_debris, start_T_sequence[i] + 5.0, t_f);
			double opt_index = N * dv;
			//double opt_index =  dv;
			if (opt_index < opt_index_min)
			{
				opt_index_min = opt_index;
				dv_min = dv;
				T_now = start_T_sequence[i];
				Ts = start_T_sequence[i] + 5.0;
				Tf = t_f;
			}
		}
		//if (dv_min > 1000.0) continue;

		if (Tf > end_epoch) continue;

		DP_info_struct dp_info_temp;
		dp_info_temp.T_single_mission = Tf;
		dp_info_temp.dv_single_mission = dv_min;
		dp_info_temp.opt_index = opt_index_min;
		dp_info[0].emplace_back(dp_info_temp);
	}


	//from NO.2 start
	for (int layer = 1; layer < sequence.size() - 1; layer++)
	{
		int now_debris = sequence[layer];
		int target_debris = sequence[layer + 1];

		vector<double> init_time;
		init_time.reserve(50);
		int size_dp_info = dp_info[layer - 1].size();
		for (int i = 0; i < size_dp_info; i++)
		{
			for (int j = 0; j < tf_num; j++)
			{
				double time_temp = dp_info[layer - 1][i].T_single_mission + 6.0 + j * t_step_tf;
				//if (time_temp > end_epoch || dp_info[layer - 1][i].dv_single_mission > 1000.0) continue;
				if (time_temp > end_epoch) continue;
				auto if_in_vector = find(init_time.begin(), init_time.end(), time_temp);
				if (if_in_vector == init_time.end())
					init_time.push_back(time_temp);
			}
		}
		//int n = init_time.size();
		for (auto t_f : init_time)
		{
			//double t_f = init_time[i];
			double opt_index_min = 1.0e8;
			double T_now, Ts, Tf;
			double dv_min = 1.0e6;
			int iter_min;
			for (int i = 0; i < size_dp_info; i++)
			{
				double t_now = dp_info[layer - 1][i].T_single_mission;
				double t_s = t_now + 5.0;
				double t_transfer = t_f - t_now - 6.0;
				if (!(0.0 - 1.0e-6 < t_transfer && t_transfer < 24.0 + 1.0e-6)) continue;
				double t_temp = t_transfer / t_step_tf;
				t_temp += 1.0e-4;
				t_temp = t_temp - int(t_temp);
				if ((t_temp) > 1.0e-3) continue;
				double dv = Dv_ij(now_debris, target_debris, t_s, t_f);
				double opt_index = (N - layer) * dv + dp_info[layer - 1][i].opt_index;
				//double opt_index = dv + dp_info[layer - 1][i].opt_index;
				if (opt_index < opt_index_min)
				{
					opt_index_min = opt_index;
					dv_min = dv;
					T_now = t_now;
					Ts = t_s;
					Tf = t_f;
					iter_min = i;
				}
			}

			//if(dv_min > 1000.0) continue;

			DP_info_struct dp_info_temp;
			dp_info_temp.last_iter = iter_min;
			dp_info_temp.T_single_mission = Tf;
			dp_info_temp.dv_single_mission = dv_min;
			dp_info_temp.opt_index = opt_index_min;
			dp_info[layer].emplace_back(dp_info_temp);
		}
	}

	//select minimum m0
	//vector<DP_info_struct_single_mission> dp_info_single_mission;
	//dp_info_single_mission.reserve(50);
	double min_m0 = 1.0e10;
	int counter_min=0;
	for (int counter_final1 = 0; counter_final1 < dp_info[N - 1].size(); counter_final1++)
	{
		if(dp_info[N - 1][counter_final1].opt_index < min_m0)
		{
			min_m0 = dp_info[N - 1][counter_final1].opt_index;
			counter_min = counter_final1;
		}
	}
	int counter_final1 = counter_min;
		DP_info_struct_single_mission dp_info_temp;
		vector < double > T_reverse, Dv_reverse;
		T_reverse.reserve(40); Dv_reverse.reserve(20);
		int counter_final = counter_final1;
		for (int i = N - 1; i > -1; i--)
		{
			auto T_temp = dp_info[i][counter_final].T_single_mission;
			auto dv_temp = dp_info[i][counter_final].dv_single_mission;
			if (i != N - 1) T_reverse.push_back(T_temp + 5.0);
			T_reverse.push_back(T_temp);
			Dv_reverse.push_back(dv_temp);
			counter_final = dp_info[i][counter_final].last_iter;
		}
		T_reverse.push_back(start_epoch + 5.0);
		T_reverse.push_back(start_epoch);

		reverse(T_reverse.begin(), T_reverse.end());
		reverse(Dv_reverse.begin(), Dv_reverse.end());

		dp_info_temp.T_single_mission.swap(T_reverse);
		mp_calc(Dv_reverse, dp_info_temp.m0);
		dp_info_temp.dv_single_mission.swap(Dv_reverse);
		//if(dp_info_temp.m0 < 7000)
		//{
		//	dp_info_single_mission.emplace_back(dp_info_temp);
		//}
	



	return dp_info_temp;
}

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
vector<DP_info_struct_single_mission> DP_optimization_single_mission(const vector<int>& sequence, double start_epoch, double end_epoch, std::vector<double>& T_single_mission)
{
	double step_init = 8.0; //days
	int num_step = 3;
	double t_step_tf = 8.0;
	int tf_num = 0;           //discrete num
	double loss_factor = 0.3;
	for (double tf_transfer_temp = 6.0; tf_transfer_temp < 30.0 + 1.0e-1; tf_transfer_temp += t_step_tf) tf_num++;
	int N = sequence.size() - 1;  //decision num

	vector<vector<DP_info_struct>> dp_info(N); //optimum_state
	for (int i = 0; i < N; i++) dp_info[i].reserve(50);

	vector<double> start_T_sequence;
	for (int i = 0; i < num_step; i++)
	{
		start_T_sequence.push_back(start_epoch + i * step_init);
	}

	if(sequence.size()==1)
	{
		vector<DP_info_struct_single_mission> temp;
		temp.resize(1);
		temp[0].m0 = 2000;
		temp[0].T_single_mission.push_back(start_epoch);
		temp[0].dv_single_mission.push_back(0.0);
		return temp;
	}
	
	//each end epoch
	int now_debris = sequence[0];
	int target_debris = sequence[1];
	int counter = 0;

	//determine x0
	set<double> init_time;
	for (int i = 0; i < start_T_sequence.size(); i++)
	{
		for (int j = 0; j < tf_num; j++)
		{
			init_time.insert(start_T_sequence[i] + 6.0 + j * t_step_tf);
		}
	}
	for (auto iter = init_time.begin(); iter != init_time.end(); ++iter)
	{
		double t_f = *iter;
		double opt_index_min = 1.0e6;
		double T_now, Ts, Tf;
		double dv_min = 1.0e6;
		for (int i = 0; i < start_T_sequence.size(); i++)
		{
			double t_transfer = t_f - start_T_sequence[i] - 6.0;
			if (!(0.0 - 1.0e-6 < t_transfer && t_transfer < 24.0 + 1.0e-6)) continue;
			double t_temp = t_transfer / t_step_tf;
			t_temp = t_temp - round(t_temp);
			if (fabs(t_temp) > 1.0e-3) continue;
			double dv = Dv_ij(now_debris, target_debris, start_T_sequence[i] + 5.0, t_f);
			double opt_index = N * dv;

			if (opt_index < opt_index_min)
			{
				opt_index_min = opt_index;
				dv_min = dv;
				T_now = start_T_sequence[i];
				Ts = start_T_sequence[i] + 5.0;
				Tf = t_f;
			}
		}
		//if (dv_min > 1000.0) continue;

		DP_info_struct dp_info_temp;
		dp_info_temp.T_single_mission=Tf;
		dp_info_temp.dv_single_mission = dv_min;
		dp_info_temp.opt_index = opt_index_min;
		dp_info[0].emplace_back(dp_info_temp);
	}


	//from NO.2 start
	for (int layer = 1; layer < sequence.size() - 1; layer++)
	{
		int now_debris = sequence[layer];
		int target_debris = sequence[layer + 1];

		vector<double> init_time;
		init_time.reserve(50);
		int size_dp_info = dp_info[layer - 1].size();
		for (int i = 0; i < size_dp_info; i++)
		{
			for (int j = 0; j < tf_num; j++)
			{
				double time_temp = dp_info[layer - 1][i].T_single_mission + 6.0 + j * t_step_tf;
				//if (time_temp > end_epoch || dp_info[layer - 1][i].dv_single_mission > 1000.0) continue;
				if (time_temp > end_epoch ) continue;
				auto if_in_vector = find(init_time.begin(), init_time.end(), time_temp);
				if(if_in_vector == init_time.end())
				init_time.push_back(time_temp);
			}
		}
		//int n = init_time.size();
		for (auto t_f: init_time)
		{
			//double t_f = init_time[i];
			double opt_index_min = 1.0e8;
			double T_now, Ts, Tf;
			double dv_min = 1.0e6;
			int iter_min;
			for (int i = 0; i < size_dp_info; i++)
			{
				double t_now = dp_info[layer - 1][i].T_single_mission;
				double t_s = t_now + 5.0;
				double t_transfer = t_f - t_now - 6.0;
				if (!(0.0 - 1.0e-6 < t_transfer && t_transfer < 24.0 + 1.0e-6)) continue;
				double t_temp = t_transfer / t_step_tf;
				t_temp += 1.0e-4;
				t_temp = t_temp - int(t_temp);
				if ((t_temp) > 1.0e-3) continue;
				double dv = Dv_ij(now_debris, target_debris, t_s, t_f);
				double opt_index = (N - layer) * dv + dp_info[layer - 1][i].opt_index;

				if (opt_index < opt_index_min)
				{
					opt_index_min = opt_index;
					dv_min = dv;
					T_now = t_now;
					Ts = t_s;
					Tf = t_f;
					iter_min = i;
				}
			}

			//if(dv_min > 1000.0) continue;

			DP_info_struct dp_info_temp;
			dp_info_temp.last_iter = iter_min;
			dp_info_temp.T_single_mission=Tf;
			dp_info_temp.dv_single_mission = dv_min;
			dp_info_temp.opt_index = opt_index_min;
			dp_info[layer].emplace_back(dp_info_temp);
		}
	}

	//select minimum m0
	vector<DP_info_struct_single_mission> dp_info_single_mission;
	dp_info_single_mission.reserve(50);
	for (int counter_final1 = 0; counter_final1 < dp_info[N - 1].size(); counter_final1++)
	{
		DP_info_struct_single_mission dp_info_temp;
		vector < double > T_reverse, Dv_reverse;
		T_reverse.reserve(40); Dv_reverse.reserve(20);
		int counter_final = counter_final1;
		for (int i = N - 1; i > -1; i--)
		{
			auto T_temp = dp_info[i][counter_final].T_single_mission;
			auto dv_temp = dp_info[i][counter_final].dv_single_mission;
			if(i!=N-1) T_reverse.push_back(T_temp+5.0);
			T_reverse.push_back(T_temp);
			Dv_reverse.push_back(dv_temp);
			counter_final = dp_info[i][counter_final].last_iter;
		}
		T_reverse.push_back(start_epoch + 5.0);
		T_reverse.push_back(start_epoch );

		reverse(T_reverse.begin(), T_reverse.end());
		reverse(Dv_reverse.begin(), Dv_reverse.end());

		dp_info_temp.T_single_mission.swap(T_reverse);
		mp_calc(Dv_reverse, dp_info_temp.m0);
		dp_info_temp.dv_single_mission.swap(Dv_reverse);
		if(dp_info_temp.m0 < 8000)
		{
			dp_info_single_mission.emplace_back(dp_info_temp);
		}
	}

	return dp_info_single_mission;
}

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
double DP_optimization_all_mission(const vector<vector<int>>& sequence, 
	                               std::vector< std::vector<double>>& T_all,
	                               std::vector< std::vector<double>>& dv_all, 
	                               double opt_min_top)
{
	//  input parameter
	const double step_init = 6.0;   // mission tf step days
	const int num_step_init = 10;   // mission tf step number
	const double t_step_tf = 6.0;   // tf step days

	// indirect parameter (do not change)
	int tf_num = 0;           //tf num
	for (double tf_transfer_temp = 6.0; tf_transfer_temp < 30.0 + 1.0e-1; tf_transfer_temp += t_step_tf) tf_num++;
	int mission_number = sequence.size();   // mission number

	double start_epoch = 0.0;
	double end_epoch = 2947.0;

	//state space 0:mission 1:debris 2:end epoch (from late to early) 3:current epoch (from late to early) ;4-d
	vector<vector<vector<vector<DP_info_struct>>>> dp_info(mission_number); 
	for(int mission=0; mission <mission_number;mission++)
	{
		int single_mission_number = sequence[mission].size();
		dp_info[mission].resize(single_mission_number);
		//for (int i = 0; i < single_mission_number; i++) dp_info[mission][i].reserve(30);
	}

	
	// start dp: determin end epoch  
	int end_single_mission_number = sequence[mission_number-1].size();
	dp_info[mission_number - 1][end_single_mission_number - 1].resize(num_step_init);
	for (int i = 0; i < num_step_init; i++)
	{
		double Tf = end_epoch - i * step_init;
		DP_info_struct dp_info_temp;
		dp_info_temp.last_iter = -1;
		dp_info_temp.T_single_mission = Tf;
		dp_info_temp.dv_single_mission = 0.0;
		dp_info_temp.opt_index = 0.0;
		dp_info_temp.end_epcoh = Tf;
		dp_info_temp.m0 = 2030.0;
		dp_info[mission_number-1][end_single_mission_number-1][i].emplace_back(dp_info_temp);
	}

	//from NO.1 start
	for (int mission_id = mission_number - 1; mission_id > -1; mission_id--)
	{
		for (int debris_id = sequence[mission_id].size() - 1; debris_id > -1; debris_id--)
		{
			//layer--;  //start from 1
			if (mission_id == 0 && debris_id == 0) break;

			if(dp_info[mission_id][debris_id].size() == 0) break;

			int counter_size = 0;
			for(int i = 0;i< dp_info[mission_id][debris_id].size();i++)
			{
				counter_size += dp_info[mission_id][debris_id][i].size();
			}
			if(counter_size==0)break;
			
			// now_debris target_debris
			int now_debris = -1;
			int target_debris = sequence[mission_id][debris_id];
			if (debris_id == 0)
			{
				now_debris = sequence[mission_id - 1].back();
			}
			else
			{
				now_debris = sequence[mission_id][debris_id - 1];
			}

			//next layer all epoch
			vector<double> next_time_all;   next_time_all.reserve(300);

			if(debris_id == 0)
			{
				//new mission;
				double t_min = 1.0e10, t_max=-1.0;
				auto& dp_info_temp = dp_info[mission_id][debris_id];
				for(auto &a: dp_info_temp)
				{
					for (auto& b : a)
					{
						double t = b.T_single_mission;
						if (t < t_min)  t_min = t;
						if (t > t_max)  t_max = t;
					}
				}
				t_min = t_min - 36.0 - (num_step_init - 1) * step_init;
				if (t_min < 1.0) t_min = 1.0;
				t_max = t_max - 36.0;
				//t_min = dp_info[mission_id][debris_id].back().back().T_single_mission - 36.0 - (num_step_init - 1) * step_init;
				//t_max = dp_info[mission_id][debris_id][0][0].T_single_mission - 36.0;
				
				for(double t_temp = t_max; t_temp > t_min - 1.0e-4; t_temp -= step_init)
				{
					next_time_all.push_back(t_temp);
				}
				int single_mission_number = sequence[mission_id - 1].size();
				dp_info[mission_id - 1][single_mission_number - 1].resize(next_time_all.size());

				for(int i = 0 ;i< next_time_all.size();i++)
				{
					double t_temp = next_time_all[i];
					DP_info_struct dp_info_temp;
					dp_info_temp.last_iter = -1;
					dp_info_temp.T_single_mission = t_temp;
					dp_info_temp.dv_single_mission = 0.0;
					dp_info_temp.opt_index = 0.0;
					dp_info_temp.end_epcoh = t_temp;
					dp_info_temp.m0 = 2030.0;
					int single_mission_number = sequence[mission_id - 1].size();
					dp_info[mission_id - 1][single_mission_number - 1][i].emplace_back(dp_info_temp);
				}
			}
			else
			{
				//same mission
				double t_min = 1.0e10, t_max = -1.0;
				auto& dp_info_temp = dp_info[mission_id][debris_id];
				for (auto& a : dp_info_temp)
				{
					for (auto& b : a)
					{
						double t = b.T_single_mission;
						if (t < t_min)  t_min = t;
						if (t > t_max)  t_max = t;
					}
				}
				t_min = t_min - 6.0 - (tf_num - 1) * t_step_tf;
				t_max = t_max - 6.0;
				if (t_min < 1.0) t_min = 1.0;
				//t_min = dp_info[mission_id][debris_id].back().back().T_single_mission - 6.0 - (tf_num - 1) * t_step_tf;
				//t_max = dp_info[mission_id][debris_id][0][0].T_single_mission - 6.0;

				for (double t_temp = t_max; t_temp > t_min - 1.0e-4; t_temp -= t_step_tf)
				{
					next_time_all.push_back(t_temp);
				}
				dp_info[mission_id][debris_id-1].resize(dp_info[mission_id][debris_id].size());
				auto& a  = dp_info[mission_id][debris_id - 1];
				for(int x =0; x<a.size();x++)
				{
					a[x].reserve(50);
				}
			}

			vector<vector<double>> dv_database; //store what need to use in calculate dv
			if(debris_id > 0)
			{
				//same mission
				vector<double> T_departure = next_time_all; reverse(T_departure.begin(), T_departure.end());
				dv_database.resize(T_departure.size());
				for(int i =0; i< dv_database.size();i++)
				{
					dv_database[i].resize(tf_num);
					double t_now = T_departure[i];
					for(int j =0; j< tf_num;j++)
					{
						dv_database[i][j] = Dv_ij(now_debris, target_debris, t_now + 5.0, t_now + 6.0 + j * t_step_tf);
					}
				}

				vector<double> end_epoch_list(dp_info[mission_id][debris_id].size());
				for(int i =0; i< end_epoch_list.size();i++)
				{
					if (dp_info[mission_id][debris_id][i].size() == 0) end_epoch_list[i] = -1.0;
					else
					end_epoch_list[i] = dp_info[mission_id][debris_id][i][0].end_epcoh;
				}


				for(int end_iter = 0; end_iter < end_epoch_list.size(); end_iter++)
				{
					double end_epoch_now = end_epoch_list[end_iter];
					if(end_epoch_now < 0.0) continue;
					int size_dp_info = dp_info[mission_id][debris_id][end_iter].size();
					//for all t_f, select best optimization index
					for (auto t_f_minor_1 : next_time_all)
					{
						double opt_index_min = opt_min_top;
						double min_end_epoch, m0;
						double Tf;
						double dv_min = 5.0e2;
						int iter_min = -1;
						for (int i = 0; i < size_dp_info; i++)
						{
							double t_now = dp_info[mission_id][debris_id][end_iter][i].T_single_mission;
							double t_s = t_f_minor_1 + 5.0;
							double t_transfer = t_now - (t_f_minor_1 + 6.0);
							if (!(0.0 - 1.0e-6 < t_transfer && t_transfer < 24.0 + 1.0e-6)) continue;
							double t_temp = t_transfer / t_step_tf;
							t_temp += 1.0e-4;
							t_temp = t_temp - int(t_temp);
							if ((t_temp) > 1.0e-3) continue;

							int row = (int) (t_f_minor_1 - next_time_all.back())/ t_step_tf;
							int colloum = (int) t_transfer / t_step_tf;
							double dv = dv_database[row][colloum];
							double m0_new = (dp_info[mission_id][debris_id][end_iter][i].m0) * exp(dv / 340.0 / 9.80665) + 30.0;
							double opt_index = (m0_new - 2000.0) * (m0_new - 2000.0) * 2.0e-6;
							if (opt_index < opt_index_min)
							{
								opt_index_min = opt_index;
								min_end_epoch = dp_info[mission_id][debris_id][end_iter][i].end_epcoh;
								m0 = m0_new;
								dv_min = dv;
								Tf = t_f_minor_1;
								iter_min = i;
							}
						}

						if (iter_min > -1)
						{
							DP_info_struct dp_info_temp;
							dp_info_temp.last_iter = iter_min;
							dp_info_temp.T_single_mission = Tf;
							dp_info_temp.dv_single_mission = dv_min;
							dp_info_temp.opt_index = opt_index_min;
							dp_info_temp.end_epcoh = min_end_epoch;
							dp_info_temp.m0 = m0;
							dp_info[mission_id][debris_id - 1][end_iter].emplace_back(dp_info_temp);
						}
						//else
						//{
						//	std::cout << "Wrong!" << std::endl;
						//}
					}
				}
			}
		}

	}

	if (dp_info[0][0].size() == 0)
	{
		T_all.clear();
		return 1.0e10;
	}
	int counter_size = 0;
	for (int i = 0; i < dp_info[0][0].size(); i++)
	{
		counter_size += dp_info[0][0][i].size();
	}
	if (counter_size == 0)
	{
		T_all.clear();
		return 1.0e10;
	}

	
	/*DP for different missions*/
	vector<vector<DP_info_mission>> missions_dp_info(mission_number);
	for (int mission = 0; mission < mission_number; mission++)
	{
		int temp= dp_info[mission][0].size();
		int mission_counter = 0;
		for(int i =0; i < temp; i++) mission_counter += dp_info[mission][0][i].size();
		missions_dp_info[mission].resize(mission_counter);
	}
	int counter = 0;
	for(int i = 0; i < dp_info[mission_number-1][0].size(); i++)
	{
		for(int j = 0; j < dp_info[mission_number - 1][0][i].size(); j++)
		{
			missions_dp_info[mission_number - 1][counter].epcoh = dp_info[mission_number - 1][0][i][j].T_single_mission;
			missions_dp_info[mission_number - 1][counter].iter = -1;
			missions_dp_info[mission_number - 1][counter].end_epoch_iter = i;
			missions_dp_info[mission_number - 1][counter].current_epoch_iter = j;
			missions_dp_info[mission_number - 1][counter].opt_index = dp_info[mission_number - 1][0][i][j].opt_index;
			counter++;
		}
	}

	
	for(int mission = mission_number-2; mission > -1; mission--)
	{
		int counter = 0;
		for (int iter_epoch = 0; iter_epoch < dp_info[mission][0].size(); iter_epoch++)
		{
			for (int i = 0; i < dp_info[mission][0][iter_epoch].size(); i++)
			{
				missions_dp_info[mission][counter].epcoh = dp_info[mission][0][iter_epoch][i].T_single_mission;
				int iter_min = -1;
				double opt_min = 1.0e20;
				for (int j = 0; j < missions_dp_info[mission + 1].size(); j++)
				{
					double t_this_end = dp_info[mission][0][iter_epoch][i].end_epcoh;
					double t_next_start = missions_dp_info[mission + 1][j].epcoh;
					if (t_next_start - t_this_end < 35.0) continue; //mission time constraint
					double opt_temp = dp_info[mission][0][iter_epoch][i].opt_index + missions_dp_info[mission + 1][j].opt_index;
					if (opt_temp < opt_min)
					{
						opt_min = opt_temp;
						iter_min = j;
					}
				}
				missions_dp_info[mission][counter].iter = iter_min;
				missions_dp_info[mission][counter].end_epoch_iter = iter_epoch;
				missions_dp_info[mission][counter].current_epoch_iter = i;
				missions_dp_info[mission][counter].opt_index = opt_min;
				counter++;
			}

		}

		
	}

	//select minimum m0
	double opt_min_each_mission = 1.0e20;
	int counter_min = -1;
	for (int counter_final1 = 0; counter_final1 < missions_dp_info[0].size(); counter_final1++)
	{
		if (missions_dp_info[0][counter_final1].opt_index < opt_min_each_mission)
		{
			opt_min_each_mission = missions_dp_info[0][counter_final1].opt_index;
			counter_min = counter_final1;
		}
	}
	vector<int> min_iter_each_mission; //each mission minimum opt iter
	int mission_id = 0;
	while (true)
	{
		min_iter_each_mission.push_back(counter_min);
		counter_min = missions_dp_info[mission_id][counter_min].iter;
		mission_id++;
		if (counter_min < 0) break;
	}


	//obtain dv t
	dv_all.resize(mission_number); vector<double> dv_single;
	T_all.resize(mission_number); vector<double> T_single;
	for(int mission = 0; mission < mission_number; mission++)
	{
		int end_epoch_iter = missions_dp_info[mission][min_iter_each_mission[mission]].end_epoch_iter;
		int current_epoch_iter = missions_dp_info[mission][min_iter_each_mission[mission]].current_epoch_iter;
		auto dp_info_now = &dp_info[mission][0][end_epoch_iter][current_epoch_iter];
		
		int debris_layer_temp = 0;
		while (true)
		{
			debris_layer_temp++;
			double dv_temp = dp_info_now->dv_single_mission;
			double t_temp = dp_info_now->T_single_mission;
			if (dv_temp < 1.0e-10) 
			{
				T_single.emplace_back(t_temp); //last one
			}
			else
			{
				T_single.emplace_back(t_temp);
				T_single.emplace_back(t_temp + 5.0);
				dv_single.emplace_back(dv_temp);
			}

			if (dp_info_now->last_iter < 0) break;

			dp_info_now = &dp_info[mission][debris_layer_temp][end_epoch_iter][dp_info_now->last_iter];
		}
		dv_all[mission].swap(dv_single);
		dv_single.clear();
		T_all[mission].swap(T_single);
		T_single.clear();
	}

	return opt_min_each_mission + 55.0 * mission_number;
}