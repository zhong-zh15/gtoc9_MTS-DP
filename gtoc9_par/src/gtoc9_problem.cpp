#include "gtoc9_problem.h"

#include <iostream>


#include "model.h"
//#include "nlopt.hpp"
//#include "PSO.h"

const double p_factor = 1e2;
const int num_test = 0;
const double loss_factor = 1.00;


/****************************************************************************
* 函数名   : load_input()
* 功  能   : 根据几次速度增量 计算mp，和初始质量
****************************************************************************/
double mp_calc(const std::vector<double>& dv, double& m0)
{
	double mf = 2000; //kg
	double m_temp = mf;
	double mp = 0.0;
	for (int i = dv.size() - 1; i > -1; i--)
	{
		m_temp += 30;
		double mp_temp = m_temp * exp(dv[i] / 340.0 / 9.80665) - m_temp;
		m_temp = m_temp + mp_temp;
		mp += mp_temp;
	}
	m_temp += 30;
	m0 = m_temp;
	return mp;
}


vector<double> divide_missions(const vector<vector<int>>& debris_sequence)
{
	int total_num = 0;
	vector<int> each_mission(debris_sequence.size());
	for (int i = 0; i < debris_sequence.size(); i++)
	{
		each_mission[i] = debris_sequence[i].size();
		total_num += each_mission[i];
	}

	double step_time = 2947.0 / total_num;


	//vector<double> endtime_each_misssion(debris_sequence.size());

	//for (int i = 0; i < debris_sequence.size(); i++) { endtime_each_misssion[i] = step_time * each_mission[i]; }

	//for (int i = 1; i < debris_sequence.size(); i++)
	//{
	//	endtime_each_misssion[i] += endtime_each_misssion[i - 1];
	//}

	//for (int i = 0; i < debris_sequence.size() - 1; i++)
	//{
	//	each_mission[i] = each_mission[i] - 35.0;
	//}


	vector<double> endtime_each_misssion(debris_sequence.size());
	for(int i =0; i< debris_sequence.size(); i++)
	{
		endtime_each_misssion[i] = (2947.0 + 35.0) / debris_sequence.size() * (i+1)-35.0;
	}
	endtime_each_misssion.back() = 2947.0;
	return endtime_each_misssion;
}



void map_input_gtoc9(vector<double>& X)
{
	//X[0] = X[0] * 500;
	for (int i = 0; i < X.size(); i++)
	{
		X[i] = X[i] * 2947.0;
		//X[i] += X[i - 1];
	}
}

void inverse_map_input_gtoc9(vector<double>& X)
{
	for (int i = 0; i < X.size(); i++)
	{
		X[i] = X[i] / 2947.0;
	}
}

double penalty_gtoc9(const vector<vector<double>>& T)
{
	double penalty = 0.0;
	double t_w = 5.0;   //再次出发最小停留时间
	double t_next_mission = 35.0; //下次到达最大任务时间
	double tf_high = 30.0; //
	double tf_ts_low = 1.0;  //出发到达的最小间距时间
	double tnow = 0.0;

	for (int mission = 0; mission < T.size(); mission++)
	{
		const vector<double>& X = T[mission];  //single mission

		//time conflict
		for (int i = 1; i < X.size(); i += 2)
		{
			tnow = X[i - 1];                                               //上一个碎片的到达时刻

			double penalt_temp = (X[i + 1] - tnow) - tf_high;              //满足下次到达最大任务时间 tf_i+1 <tf_i + 30 
			if (penalt_temp > 0.0)
				penalty += p_factor * penalt_temp;

			double penalt_temp2 = t_w + tnow - X[i];                       //满足再次出发最小停留时间 tf_i + 5 <ts_i+1
			if (penalt_temp2 > 0.0)
				penalty += p_factor * penalt_temp2;

			double penalt_temp3 = X[i] + tf_ts_low - X[i + 1];             //满足出发到达的最小间距时间 ts_i+1 + 1 <tf_i+1
			if (penalt_temp3 > 0.0)
				penalty += p_factor * penalt_temp3;
		}
	}
	//next mission
	for (int mission = 1; mission < T.size(); mission++)
	{
		double penalt_temp = T[mission - 1].back() + t_next_mission - T[mission][0];
		if (penalt_temp > 0.0) penalty += p_factor * penalt_temp;
	}

	return penalty;
}




/****************************************************************************
* Function     : Obj_gtoc9()
* Discription  :
*                input: time squence(0~1 -> 0~2947 days)
*                para:  debris squence (Opt_info_gtoc9)
*                ouput: cost function of GTOC 9
****************************************************************************/
double Obj_gtoc9(const std::vector<double>& X, std::vector<double>& grad, void* f_data)
{
	vector<double> T = X;
	map_input_gtoc9(T);
	Opt_info_gtoc9* para = (Opt_info_gtoc9*)f_data; //transfer to Opt_info_gtoc9
	vector<vector<int>>& debris_squence = para->debris_squence;

	int mission_number = para->debris_squence.size();

	//T_all, dv for all missions
	std::vector<vector<double>> T_all(mission_number), dv_missions(mission_number);

	int size = 0; for (int i = 0; i < mission_number; i++) size += debris_squence[i].size();
	if (T.size() == size)
	{
		//if time squence only contains tf (i.e. sum_mission (n) )
		auto iter = T.begin();
		for (int i = 0; i < mission_number; i++)
		{
			T_all[i].assign(iter, iter + debris_squence[i].size());
			for (int j = T_all[i].size() - 1; j > 0; j--) T_all[i].insert(T_all[i].begin() + j, *(T_all[i].begin() + (j - 1)));
			iter += debris_squence[i].size();
			Dv_All(T_all[i], debris_squence[i], dv_missions[i]);
		}
	}
	else if (T.size() == 2 * size - mission_number)
	{
		// if time squence contains ts tf  (i.e. sum_mission (2n-1) )
		auto iter = T.begin();
		for (int i = 0; i < mission_number; i++)
		{
			T_all[i].assign(iter, iter + (2 * debris_squence[i].size() - 1));
			iter += (2 * debris_squence[i].size() - 1);
			Dv_All(T_all[i], debris_squence[i], dv_missions[i]);
		}
	}
	else
	{
		std::cout << "time squence number wrong!";
	}

	//Calculate cost function
	double total_socre = 0.0;
	for (int mission_id = 0; mission_id < dv_missions.size(); mission_id++)
	{
		//Calculate mass
		double m0;
		double mp = mp_calc(dv_missions[mission_id], m0);
		total_socre += 55.0 + 2.0e-6 * (m0 - 2000.0) * (m0 - 2000.0);
	}

	double penalty = penalty_gtoc9(T_all);

	return total_socre + penalty;
}


vector<double> out_gtoc9_index_m0(const std::vector<double>& X, void* f_data)
{
	vector<double> T = X;
	//map_input_gtoc9(T);
	Opt_info_gtoc9* para = (Opt_info_gtoc9*)f_data; //transfer to Opt_info_gtoc9
	vector<vector<int>>& debris_squence = para->debris_squence;

	int mission_number = para->debris_squence.size();

	//T_all, dv for all missions
	std::vector<vector<double>> T_all(mission_number), dv_missions(mission_number);

	int size = 0; for (int i = 0; i < mission_number; i++) size += debris_squence[i].size();
	if (T.size() == size)
	{
		//if time squence only contains tf (i.e. sum_mission (n) )
		auto iter = T.begin();
		for (int i = 0; i < mission_number; i++)
		{
			T_all[i].assign(iter, iter + debris_squence[i].size());
			for (int j = T_all[i].size() - 1; j > 0; j--) T_all[i].insert(T_all[i].begin() + j, *(T_all[i].begin() + (j - 1)));
			iter += debris_squence[i].size();
			Dv_All(T_all[i], debris_squence[i], dv_missions[i]);
		}
	}
	else if (T.size() == 2 * size - mission_number)
	{
		// if time squence contains ts tf  (i.e. sum_mission (2n-1) )
		auto iter = T.begin();
		for (int i = 0; i < mission_number; i++)
		{
			T_all[i].assign(iter, iter + (2 * debris_squence[i].size() - 1));
			iter += (2 * debris_squence[i].size() - 1);
			Dv_All(T_all[i], debris_squence[i], dv_missions[i]);
		}
	}
	else
	{
		std::cout << "time squence number wrong!";
	}

	//Calculate cost function
	double total_socre = 0.0;
	for (int mission_id = 0; mission_id < dv_missions.size(); mission_id++)
	{
		//Calculate mass
		double m0;
		double mp = mp_calc(dv_missions[mission_id], m0);
		total_socre += 55.0 + 2.0e-6 * (m0 - 2000.0) * (m0 - 2000.0);
	}

	double penalty = penalty_gtoc9(T_all);

	vector<double> m0_vector;

	for(int i = 0; i< dv_missions.size(); i++)
	{
		//Calculate mass
		double m0;
		double mp = mp_calc(dv_missions[i], m0);
		m0_vector.push_back(m0);
		cout << "Mission ID： " << i+1 << "  m0 : " << m0 << endl;
	}

	cout << "score :" << total_socre << endl;
	return m0_vector;
}


/****************************************************************************
* Function     : Obj_gtoc9_single_mission_single_debris()
* Discription  :
*                input: time squence(0~1 -> 0~2947 days)
*                para:  debris squence (Opt_info_gtoc9)
*                ouput: cost function of GTOC 9
****************************************************************************/
double Obj_gtoc9_single_mission_single_debris(const std::vector<double>& X, std::vector<double>& grad, void* f_data)
{
	vector<double> T_temp = X;
	map_input_gtoc9(T_temp);
	Opt_info_gtoc9* para = (Opt_info_gtoc9*)f_data; //transfer to Opt_info_gtoc9
	vector<vector<int>>& debris_squence = para->debris_squence;
	vector<vector<double>>& time_sequence = para->time_sequence;
	int single_mission = para->mission;
	double end_epoch = para->end_peoch;
	int N = para->last_debris_num;

	//generate T sequence
	vector<double> T;
	for(int mission = 0; mission < debris_squence.size();mission++)
	{
		T.insert(T.end(), time_sequence[mission].begin(), time_sequence[mission].end());
		if(mission == single_mission) T.insert(T.end(), T_temp.begin(), T_temp.end());
	}

	int mission_number = para->debris_squence.size();

	//T_all, dv for all missions
	std::vector<vector<double>> T_all(mission_number), dv_missions(mission_number);

	int size = 0; for (int i = 0; i < mission_number; i++) size += debris_squence[i].size();
	if (T.size() == size)
	{
		//if time squence only contains tf (i.e. sum_mission (n) )
		auto iter = T.begin();
		for (int i = 0; i < mission_number; i++)
		{
			T_all[i].assign(iter, iter + debris_squence[i].size());
			for (int j = T_all[i].size() - 1; j > 0; j--) T_all[i].insert(T_all[i].begin() + j, *(T_all[i].begin() + (j - 1)));
			iter += debris_squence[i].size();
			Dv_All(T_all[i], debris_squence[i], dv_missions[i]);
		}
	}
	else if (T.size() == 2 * size - mission_number)
	{
		// if time squence contains ts tf  (i.e. sum_mission (2n-1) )
		auto iter = T.begin();
		for (int i = 0; i < mission_number; i++)
		{
			T_all[i].assign(iter, iter + (2 * debris_squence[i].size() - 1));
			iter += (2 * debris_squence[i].size() - 1);
			Dv_All(T_all[i], debris_squence[i], dv_missions[i]);
		}
	}
	else
	{
		std::cout << "time squence number wrong!";
	}

	//Calculate cost function
	double total_socre = 0.0;
	for (int mission_id = 0; mission_id < dv_missions.size(); mission_id++)
	{
		//Calculate mass
		double m0;
		double mp = mp_calc(dv_missions[mission_id], m0);
		total_socre += 55.0 + 2.0e-6 * (m0 - 2000.0) * (m0 - 2000.0);
	}

	double penalty = penalty_gtoc9(T_all);

	//single mission time_penalty
	//double T_last =  (T_all[single_mission].back() - loss_factor * (end_epoch - (N - 1) * (end_epoch - T_all[single_mission][0]) / (T_all[single_mission].size() + N - 2 + num_test)));

	double t_now = T_all[single_mission].back();
	double t_0 = T_all[single_mission][0];
	int n = debris_squence[single_mission].size() - 1;

	double T_last = (t_now - t_0)/(end_epoch - t_0) - pow(( (n + 1.0) / (N + n)), loss_factor);

	if(T_last > 0.0)
	{
		penalty += p_factor * p_factor* 30*T_last;
	}

	//ts < 6 day time conflict
	vector<double >& T_single_mission = T_all[single_mission];
	for (int i = 1; i < T_single_mission.size(); i += 2)
	{
		double tnow = T_single_mission[i - 1];                        //上一个碎片的到达时刻

		double penalt_temp2 = -(5.1 + tnow) + T_single_mission[i];                       //满足再次出发最小停留时间 tf_i + 5 <ts_i+1
		if (penalt_temp2 > 0.0)
			penalty += p_factor * penalt_temp2;
	}

	return total_socre + penalty;
}



double Obj_gtoc9_single_mission(const std::vector<double>& X, std::vector<double>& grad, void* f_data)
{
	vector<double> T_temp = X;
	map_input_gtoc9(T_temp);
	Opt_info_gtoc9* para = (Opt_info_gtoc9*)f_data; //transfer to Opt_info_gtoc9
	vector<vector<int>>& debris_squence = para->debris_squence;
	vector<vector<double>>& time_sequence = para->time_sequence;
	int single_mission = para->mission;
	double end_epoch = para->end_peoch;
	int N = para->last_debris_num;

	//generate T sequence
	vector<double> T;
	for (int mission = 0; mission < debris_squence.size(); mission++)
	{
		if (mission == single_mission) T.insert(T.end(), T_temp.begin(), T_temp.end());
		else T.insert(T.end(), time_sequence[mission].begin(), time_sequence[mission].end());
	}

	int mission_number = para->debris_squence.size();

	//T_all, dv for all missions
	std::vector<vector<double>> T_all(mission_number), dv_missions(mission_number);

	int size = 0; for (int i = 0; i < mission_number; i++) size += debris_squence[i].size();
	if (T.size() == size)
	{
		//if time squence only contains tf (i.e. sum_mission (n) )
		auto iter = T.begin();
		for (int i = 0; i < mission_number; i++)
		{
			T_all[i].assign(iter, iter + debris_squence[i].size());
			for (int j = T_all[i].size() - 1; j > 0; j--) T_all[i].insert(T_all[i].begin() + j, *(T_all[i].begin() + (j - 1)));
			iter += debris_squence[i].size();
			Dv_All(T_all[i], debris_squence[i], dv_missions[i]);
		}
	}
	else if (T.size() == 2 * size - mission_number)
	{
		// if time squence contains ts tf  (i.e. sum_mission (2n-1) )
		auto iter = T.begin();
		for (int i = 0; i < mission_number; i++)
		{
			T_all[i].assign(iter, iter + (2 * debris_squence[i].size() - 1));
			iter += (2 * debris_squence[i].size() - 1);
			Dv_All(T_all[i], debris_squence[i], dv_missions[i]);
		}
	}
	else
	{
		std::cout << "time squence number wrong!";
	}

	//Calculate cost function
	double total_socre = 0.0;
	for (int mission_id = 0; mission_id < dv_missions.size(); mission_id++)
	{
		//Calculate mass
		double m0;
		double mp = mp_calc(dv_missions[mission_id], m0);
		total_socre += 55.0 + 2.0e-6 * (m0 - 2000.0) * (m0 - 2000.0);
	}

	double penalty = penalty_gtoc9(T_all);

	//single mission time_penalty
	//double T_last =  (T_all[single_mission].back() - loss_factor * (end_epoch - (N - 1) * (end_epoch - T_all[single_mission][0]) / (T_all[single_mission].size() + N - 2 + num_test)));

	double t_now = T_all[single_mission].back();
	double t_0 = T_all[single_mission][0];
	int n = debris_squence[single_mission].size() - 1;

	double T_last = (t_now - t_0) / (end_epoch - t_0) - pow(((n + 1.0) / (N + n)), loss_factor);

	//if (n < 4) T_last = -1.0;

	if (T_last > 0.0)
	{
		penalty += p_factor * p_factor * 30 * T_last;
	}

	return total_socre + penalty;
}

vector<double> estimate_t_sequence(const vector<vector<int>>& debris_sequence)
{
	vector<vector<double>> T_allmissions(debris_sequence.size());

	//end epoch
	auto endtime_each_misssion = divide_missions(debris_sequence);

	//start epoch for each mission
	T_allmissions[0].push_back(0.0);
	for (int i = 1; i < T_allmissions.size(); i++)
	{
		T_allmissions[i].push_back(endtime_each_misssion[i - 1] + 35.0);
	}

	//double estimate_dv(int debris_now, int target, double tnow, double& Ts, double& Tf, double end_epoch, int N)

	//t_sequence
	for (int mission = 0; mission < debris_sequence.size(); mission++)
	{
		for (int i = 1; i < debris_sequence[mission].size(); i++)
		{
			double ts, tf;
			double tnow = T_allmissions[mission].back();
			int debris_now = debris_sequence[mission][i - 1];
			int debris_target = debris_sequence[mission][i];
			estimate_dv(debris_now, debris_target, tnow, ts, tf, endtime_each_misssion[mission], debris_sequence[mission].size() - i);
			T_allmissions[mission].push_back(ts);
			T_allmissions[mission].push_back(tf);
		}
	}

	vector<double> T;
	for (auto& a : T_allmissions)
	{
		for (auto& b : a)
		{
			T.push_back(b);
		}
	}

	return T;
}


//vector<double> estimate_t_sequence_random(const vector<vector<int>>& debris_sequence)
//{
//	vector<vector<double>> T_allmissions(debris_sequence.size());
//
//	//end epoch
//	auto endtime_each_misssion = divide_missions(debris_sequence);
//
//	//start epoch for each mission
//	T_allmissions[0].push_back(0.0);
//	for (int i = 1; i < T_allmissions.size(); i++)
//	{
//		T_allmissions[i].push_back(endtime_each_misssion[i - 1] + 35.0);
//	}
//
//	vector<double> start_MJD_eachmission = { 23557.18, 23851.08, 24057.47, 24637.26, 24946.47, 25262.95, 25485.20, 25712.38, 25946.06, 26267.80 };
//	for (int i = 0; i < T_allmissions.size(); i++)
//		T_allmissions[i][0] = start_MJD_eachmission[i] - 23467.0;
//	
//	//random
//	//double offset = 30.0;
//	//T_allmissions[0][0] += offset;
//	//for (int i = 0; i < T_allmissions.size(); i++)
//	//{
//
//	//	T_allmissions[i][0] += realRand(-offset, offset);
//	//}
//
//	//double estimate_dv(int debris_now, int target, double tnow, double& Ts, double& Tf, double end_epoch, int N)
//
//	//t_sequence
//	for (int mission = 0; mission < debris_sequence.size(); mission++)
//	{
//		for (int i = 1; i < debris_sequence[mission].size(); i++)
//		{
//			double ts, tf;
//			double tnow = T_allmissions[mission].back();
//			int debris_now = debris_sequence[mission][i - 1];
//			int debris_target = debris_sequence[mission][i];
//			estimate_dv_random(debris_now, debris_target, tnow, ts, tf, endtime_each_misssion[mission], debris_sequence[mission].size() - i);
//			T_allmissions[mission].push_back(ts);
//			T_allmissions[mission].push_back(tf);
//		}
//	}
//
//	vector<double> T;
//	for (auto& a : T_allmissions)
//	{
//		for (auto& b : a)
//		{
//			T.push_back(b);
//		}
//	}
//
//	return T;
//}


//void optimize_single_mission_ts_tf(vector<vector<double>>& T_sequence, Opt_info_gtoc9& opt_info, int N)
//{
//	int mission = opt_info.mission;
//	vector<vector<int>>& debris_squence = opt_info.debris_squence;
//	vector<vector<double>>& time_sequence = opt_info.time_sequence;
//	int debris_now = debris_squence[mission][debris_squence[mission].size() - 2];
//	double Ts, Tf;
//	double end_epoch;
//	if (mission == debris_squence.size() - 1)
//	{
//		end_epoch = 2947.0;
//	}
//	else
//	{
//		end_epoch = time_sequence[mission + 1][0] - 35.0;
//	}
//	opt_info.end_peoch = end_epoch;
//	opt_info.last_debris_num = N;
//	estimate_dv(debris_now, debris_squence[mission].back(), time_sequence[mission].back(), Ts, Tf, end_epoch, N);
//	//estimate_dv_random(debris_now, debris_squence[mission].back(), time_sequence[mission].back(), Ts, Tf, end_epoch, N);
//	//NLOPT optimize
//	int num_variable = 2;                                //only ts, tf
//	int ItMax = 1e5;
//	nlopt::opt opter(nlopt::LN_SBPLX, num_variable);				//局部优化  LN_SBPLX  LN_COBYLA LN_BOBYQA
//																			  //GN_DIRECT x  GD_STOGO x  GN_AGSx GN_CRS2_LM GN_ISRESx GN_ESCH
//																			  //效果来看 LN_SBPLX、 GN_CRS2_LM  比较好
//	opter.set_min_objective(Obj_gtoc9_single_mission_single_debris, &opt_info);					//指标
//	std::vector<double> lb(num_variable), rb(num_variable), dx(num_variable);
//	//opter.add_inequality_mconstraint(ObjFun, node_order, dx);
//	double fbest = 1.0e20;
//	for (int i = 0; i < num_variable; i++)
//	{
//		lb[i] = 0.0;												//下界
//		rb[i] = 1.0;												//上界
//		dx[i] = 1.0 / 86400.0 * 200.0;									//初始步长
//	}
//	//opter.set_initial_step(dx);								    //设置初始步长
//	opter.set_lower_bounds(lb);										//设置下届
//	opter.set_upper_bounds(rb);										//设置上届
//	double tol = 1e-8;
//	//opter.set_ftol_abs(tol);
//	//opter.set_force_stop(tol);
//	opter.set_maxeval(ItMax);										//优化该次数后停止
//	vector<double> T = { Ts,Tf };
//	inverse_map_input_gtoc9(T);
//	nlopt::result res = opter.optimize(T, fbest);				 //进行优化
//
//	std::cout << fbest << std::endl;
//
//	map_input_gtoc9(T);
//	// transfer to multi missions
//	T_sequence[mission].push_back(T[0]);
//	T_sequence[mission].push_back(T[1]);
//
//	opt_info.time_sequence[mission].push_back(T[0]);
//	opt_info.time_sequence[mission].push_back(T[1]);
//
//}
//

//void optimize_single_mission(vector<vector<double>>& T_sequence, Opt_info_gtoc9& opt_info, int N)
//{
//	int mission = opt_info.mission;
//	vector<vector<int>>& debris_squence = opt_info.debris_squence;
//	vector<vector<double>>& time_sequence = opt_info.time_sequence;
//
//	double end_epoch;
//	if (mission == debris_squence.size() - 1)
//	{
//		end_epoch = 2947.0;
//	}
//	else
//	{
//		end_epoch = time_sequence[mission + 1][0] - 35.0;
//	}
//	opt_info.end_peoch = end_epoch;
//	opt_info.last_debris_num = N;
//
//	//NLOPT optimize
//	int num_variable = time_sequence[mission].size();                                //only ts, tf
//	int ItMax = 1e5;
//	nlopt::opt opter(nlopt::LN_SBPLX, num_variable);				//局部优化  LN_SBPLX  LN_COBYLA LN_BOBYQA
//																			  //GN_DIRECT x  GD_STOGO x  GN_AGSx GN_CRS2_LM GN_ISRESx GN_ESCH
//																			  //效果来看 LN_SBPLX、 GN_CRS2_LM  比较好
//	opter.set_min_objective(Obj_gtoc9_single_mission, &opt_info);					//指标
//	std::vector<double> lb(num_variable), rb(num_variable), dx(num_variable);
//	//opter.add_inequality_mconstraint(ObjFun, node_order, dx);
//	double fbest = 1.0e20;
//	for (int i = 0; i < num_variable; i++)
//	{
//		lb[i] = 0.0;												//下界
//		rb[i] = 1.0;												//上界
//		dx[i] = 1.0 / 86400.0 * 200.0;									//初始步长
//	}
//	//opter.set_initial_step(dx);								    //设置初始步长
//	opter.set_lower_bounds(lb);										//设置下届
//	opter.set_upper_bounds(rb);										//设置上届
//	double tol = 1e-8;
//	//opter.set_ftol_abs(tol);
//	//opter.set_force_stop(tol);
//	opter.set_maxeval(ItMax);										//优化该次数后停止
//	vector<double> T = T_sequence[mission];
//	inverse_map_input_gtoc9(T);
//	nlopt::result res = opter.optimize(T, fbest);				 //进行优化
//
//	std::cout << fbest << std::endl;
//
//	map_input_gtoc9(T);
//	// transfer to multi missions
//	T_sequence[mission] = T;
//	T_sequence[mission] = T;
//
//	opt_info.time_sequence[mission] = T;
//	opt_info.time_sequence[mission] = T;
//
//}
//
//void single_mission(const vector<double>& mission_start_epoch, const Opt_info_gtoc9& opt_info, vector< vector<double>>& T_sequence)
//{
//	const int mission_number = opt_info.debris_squence.size();
//
//	//inital sequence, only start debris
//	Opt_info_gtoc9 opt_info_temp;
//	opt_info_temp.debris_squence.resize(mission_number);
//	for (int mission = 0; mission < mission_number; mission++) opt_info_temp.debris_squence[mission].push_back(opt_info.debris_squence[mission][0]);
//
//	//inital T_sequence, only start epoch
//	T_sequence.clear(); T_sequence.resize(mission_number);
//	for (int mission = 0; mission < mission_number; mission++) T_sequence[mission].push_back(mission_start_epoch[mission]);
//	opt_info_temp.time_sequence = T_sequence;
//
//	// add debris one by one
//	for (int mission = 0; mission < mission_number; mission++)
//	{
//		for (int target_debris_num = 1; target_debris_num < opt_info.debris_squence[mission].size(); target_debris_num++)
//		{
//			opt_info_temp.debris_squence[mission].push_back(opt_info.debris_squence[mission][target_debris_num]);
//			opt_info_temp.mission = mission;
//			//optimize_single_mission_ts_tf(T_sequence, opt_info_temp, opt_info.debris_squence[mission].size()-target_debris_num);
//			//optimize_single_mission(T_sequence, opt_info_temp, opt_info.debris_squence[mission].size() - target_debris_num);
//		}
//		//cout << " Now is to optimize whole mission！" << std::endl;
//		//optimize_single_mission(T_sequence, opt_info_temp, 1);
//
//
//		vector<double> T_single;
//		vector<int> Debris_single = opt_info.debris_squence[mission];
//		double start_epoch = mission_start_epoch[mission];
//		double end_epoch;
//		if (mission == mission_number - 1)
//			end_epoch = 2947.0;
//		else
//			end_epoch = mission_start_epoch[mission + 1] - 35.0;
//
//		//DP_optimization(Debris_single, start_epoch, end_epoch, T_single);
//		auto& T_temp_all_missions = T_sequence;
//		T_temp_all_missions[mission] = T_single;
//		//cout << " Now is to optimize DP!!!! single mission！" << std::endl;
//		//optimize_single_mission(T_temp_all_missions, opt_info_temp, 1);
//
//	}
//
//	//optimize maximum m0 
//
//	//vector<double> T_test_allmission;
//	//for (int i = 0; i < T_sequence.size(); i++) T_test_allmission.insert(T_test_allmission.end(),
//	//	T_sequence[i].begin(), T_sequence[i].end());
//
//	//inverse_map_input_gtoc9(T_test_allmission);
//	//auto m0 = out_gtoc9_index_m0(T_test_allmission, &opt_info_temp);
//
//	//auto m0_temp = m0;
//	//sort(m0_temp.rbegin(), m0_temp.rend());
//	//for(int i = 0; i< mission_number ;i++)
//	//{
//	//	auto iter = find(m0. begin(), m0.end(), m0_temp[i]);
//	//	int mission = distance(m0.begin(), iter);
//	//	opt_info_temp.mission = mission;
//	//	optimize_single_mission(T_sequence, opt_info_temp, 1);
//	//}
//
//}
//

//double time_optimization_NLOPT(vector<double>& T_opt, const vector<vector<int>>& debris_sequence)
//{
//	Opt_info_gtoc9 JPL_opt_info;
//	JPL_opt_info.debris_squence = debris_sequence;
//
//	int num_variable = T_opt.size();                                //只优化tf
//	int ItMax = 1e8;
//	auto opter = nlopt::opt(nlopt::LN_SBPLX, num_variable);
//	//GN_DIRECT x  GD_STOGO x  GN_AGSx GN_CRS2_LM GN_ISRESx GN_ESCH
//	//效果来看 LN_SBPLX、 GN_CRS2_LM  比较好
//
//	opter.set_min_objective(Obj_gtoc9, &JPL_opt_info);					//指标
//
//	std::vector<double> lb(num_variable), rb(num_variable), dx(num_variable);
//
//	//opter.add_inequality_mconstraint(ObjFun, node_order, dx);
//
//	double fbest = 1.0e20;
//	for (int i = 0; i < num_variable; i++)
//	{
//		lb[i] = 0.0;												//下界
//		rb[i] = 1.0;												//上界
//		dx[i] = 1.0 / 86400.0 * 200.0;									//初始步长
//	}
//	//opter.set_initial_step(dx);										//设置初始步长
//	opter.set_lower_bounds(lb);										//设置下届
//	opter.set_upper_bounds(rb);										//设置上届
//	double tol = 1e-8;
//	//opter.set_ftol_abs(tol);
//	//opter.set_force_stop(tol);
//
//	opter.set_maxeval(ItMax);										//优化该次数后停止
//
//	//vector<double> T = T_opt;
//	inverse_map_input_gtoc9(T_opt);
//
//	nlopt::result res = opter.optimize(T_opt, fbest);				   //进行优化
//	map_input_gtoc9(T_opt);
//
//	std::cout << fbest << std::endl;
//	return fbest;
//}