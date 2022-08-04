/****************************************************************************
* Copyright (C), 2020-2031 清华大学航天航空学院动力学与控制实验室
* 作者: 张众 zhong-zh19@mails.tsinghua.edu.cn
*            539977562@qq.com
* 文件名: main.cpp
* 内容简述：GTOC9并行算例的主程序
*
* 文件历史：
* 版本号     日期         作者       说明
* 01       2021-05-12    张众      创建文件
****************************************************************************/

#include "main.h"

#include <chrono>
#include "DP.h"
#include "gtoc9_problem.h"


//数据集debris_data
//0碎片编号 1* 2轨道半长轴(m) 3偏心率 4倾角(rad) 
//5升交点赤经(rad) 6近地点幅角(rad) 7真近点角(rad)
double debris_data[123][8];   //全局碎片信息
//double**** debris_database_dv;
const int debrisnum = 123;
const int departure_time = 591;
const int end_time = 13;
//double** debris_database_dv_p;

double domega_debris[123];
double dOmega_debris[123];
//double dOmega_debris_notreal[123];

/****************************************************************************
* 函数名   : load_input()
* 功  能   : 读取空间碎片信息,并做预处理
****************************************************************************/
void load_input()
{
	double debris_data_temp[123][8];
	std::ifstream fin("../input_data/debris_data.txt");
	for (int i = 0; i < 123; i++)
	{
		double temp;
		fin >> temp;
		for (int j = 0; j < 7; j++)
			fin >> debris_data_temp[i][j];
		int flag = 0;
		debris_data_temp[i][7] = E2f(flag, M2E(flag, debris_data_temp[i][6], debris_data_temp[i][2]), debris_data_temp[i][2]);
	}
/*
 * 0 初始时刻
 * 1 a 单位：m
 * 2 e
 * 3 i 单位：rad
 * 4 W 单位：rad
 * 5 omega 单位：rad
 * 6 M 单位：rad
 * 7 f 单位：rad
 */

	for(int i = 0; i<123; i++)
	{
		double Re = req;// 6378137.0;													//m 地球赤道半径

		double a = debris_data_temp[i][1];
		double e = debris_data_temp[i][2];
		double inc = debris_data_temp[i][3];
		double OMEGA = debris_data_temp[i][4];
		double omega = debris_data_temp[i][5];
		double M = debris_data_temp[i][6];

		
		double p = a * (1 - e * e);
		double c2 = (Re / p) * (Re / p);
		double ci = cos(inc);
		double n = sqrt(mu / (a * a * a));
		double dOmega = -1.5 * J2 * c2 * n * ci;
		double domega = 0.75 * J2 * c2 * n * (5 * ci * ci - 1);
		double dt = (23467.0 - debris_data_temp[i][0]) * Day2Second;

		OMEGA += dOmega * dt;
		omega += domega * dt;
		M += n * dt;
		int flag = 0;
		double f = E2f(flag, M2E(flag, M, e), e);
		//数据集debris_data
		//0碎片编号 1* 2轨道半长轴(m) 3偏心率 4倾角(rad) 5升交点赤经(rad) 6近地点幅角(rad) 7真近点角(rad)
		debris_data[i][2] = a;
		debris_data[i][3] = e;
		debris_data[i][4] = inc;
		debris_data[i][5] = OMEGA;
		debris_data[i][6] = omega;
		debris_data[i][7] = f;
	}
	for (int i = 0; i < 123; i++)
	{
		domega_init(i);
		dOmega_init(i);
		/*dOmega_notreal(i);*/
	}
		

	//需要将信息转化为初始时刻即MJD23467
	fin.close();

	
}


void Test_GTOC9_results()
{
	//Nan Zhang solution
	vector < int > Nan_1 = { 120, 114 ,15, 50, 55, 113, 121, 117, 79, 25, 27, 84, 20 ,118, 87 };
	vector < int > Nan_2 = { 76, 28, 66, 58, 74, 29, 52, 64, 51, 72, 61, 107 };
	vector < int > Nan_3 = { 41, 7, 77, 45, 82, 70, 88, 85, 47, 104, 75, 18, 6 };
	vector < int > Nan_4 = { 10, 43, 8, 4, 71, 110, 73, 14, 9, 95, 93, 90, 19, 69, 21 };
	vector < int > Nan_5 = { 81, 63, 112, 119, 105, 24, 46, 96, 108, 37, 35, 32 };
	vector < int > Nan_6 = { 36, 59, 98, 11, 17, 91, 62, 122, 54, 1, 89, 0, 40 };
	vector < int > Nan_7 = { 5, 39, 12, 101, 48, 2, 38, 106, 33, 68, 26 };
	vector < int > Nan_8 = { 97, 115, 22, 102, 86, 34, 100, 30, 92, 65, 31 };
	vector < int > Nan_9 = { 13, 94, 78, 103, 111, 56, 67, 3, 42, 44, 57 };
	vector < int > Nan_10 = { 109, 80, 16, 83, 53, 49, 23, 99, 60, 116 };
	vector<vector<int>> Nan_sequence;
	Nan_sequence.push_back(Nan_1);
	Nan_sequence.push_back(Nan_2);
	Nan_sequence.push_back(Nan_3);
	Nan_sequence.push_back(Nan_4);
	Nan_sequence.push_back(Nan_5);
	Nan_sequence.push_back(Nan_6);
	Nan_sequence.push_back(Nan_7);
	Nan_sequence.push_back(Nan_8);
	Nan_sequence.push_back(Nan_9);
	Nan_sequence.push_back(Nan_10);


	//JPL solution
	//Sequence
	vector < int > JPL_1 = { 23,55,79,113,25,20,27,117,121,50,95,102,38,97 };
	vector < int > JPL_2 = { 19,115,41,26,45,82,47,85,7,2,11,77 };
	vector < int > JPL_3 = { 72,107,61,10,28,3,64,66,31,90,73,87,57,35,69,65,8,43,71,4,29 };
	vector < int > JPL_4 = { 108,24,104,119,22,75,63,112,37,32,114 };
	vector < int > JPL_5 = { 84,59,98,1,40,51,36,67,62,99,54,122,76,15 };
	vector < int > JPL_6 = { 101,48,53,5,12,39,58,13,60,74 };
	vector < int > JPL_7 = { 49,9,70,93,105,46,88,118,18,91 };
	vector < int > JPL_8 = { 86,34,100,30,92,6,110,96,81 };
	vector < int > JPL_9 = { 33,68,116,106,14,52,120,80,16,94,83,89 };
	vector < int > JPL_10 = { 44,111,56,78,0,17,109,103,42,21 };
	vector<vector<int>> JPL_sequence;
	JPL_sequence.push_back(JPL_1);
	JPL_sequence.push_back(JPL_2);
	JPL_sequence.push_back(JPL_3);
	JPL_sequence.push_back(JPL_4);
	JPL_sequence.push_back(JPL_5);
	JPL_sequence.push_back(JPL_6);
	JPL_sequence.push_back(JPL_7);
	JPL_sequence.push_back(JPL_8);
	JPL_sequence.push_back(JPL_9);
	JPL_sequence.push_back(JPL_10);

	//mission start time
	vector<double> start_MJD_eachmission = { 23557.18, 23851.08, 24057.47, 24637.26, 24946.47, 25262.95, 25485.20, 25712.38, 25946.06, 26267.80 };
	for (int i = 0; i < start_MJD_eachmission.size(); i++) start_MJD_eachmission[i] = start_MJD_eachmission[i] - 23467.0;

	//Epoch
	vector < double > JPL_rendezvous_1 = { 5.00,5.00,5.04,5.01,5.01,5.03,5.00,5.00,5.00,5.03,5.03,5.04,5.04,5.00 };
	vector < double > JPL_rendezvous_2 = { 5.00,5.02,5.02,5.00,5.04,5.00,5.05,5.02,5.07,5.03,5.02,5.00 };
	vector < double > JPL_rendezvous_3 = { 5.00,5.06,5.01,5.02,5.07,5.02,5.04,5.02,5.01,5.02,5.01,5.07,5.06,5.02,5.01,5.01,5.06,5.01,5.02,5.04,5.00 };
	vector < double > JPL_rendezvous_4 = { 5.00,6.01,6.01,6.03,6.05,6.05,6.04,6.01,6.06,6.04,5.00 };
	vector < double > JPL_rendezvous_5 = { 5.00,5.02,5.07,5.04,5.01,5.01,5.02,5.06,5.06,5.02,5.06,5.01,5.07,5.00 };
	vector < double > JPL_rendezvous_6 = { 5.00,5.02,5.01,5.04,5.07,5.02,5.01,5.02,5.02,5.00 };
	vector < double > JPL_rendezvous_7 = { 5.00,5.00,5.06,5.06,5.04,5.06,5.04,5.06,5.03,5.00 };
	vector < double > JPL_rendezvous_8 = { 5.00,5.01,5.03,5.00,5.01,5.04,5.07,5.02,5.00 };
	vector < double > JPL_rendezvous_9 = { 5.00,5.51,5.53,5.53,5.53,5.55,5.54,5.53,5.54,5.55,5.52,5.00 };
	vector < double > JPL_rendezvous_10 = { 5.00,5.54,5.50,5.50,5.52,5.52,5.54,5.53,5.52,5.00 };

	vector < double > JPL_transfer_1 = { 24.86,24.98,22.42,24.99,0.29,10.63,25.00,2.70,1.51,1.41,24.67,24.31,5.86 };
	vector < double > JPL_transfer_2 = { 24.93,0.28,0.73,0.39,17.07,1.61,22.42,2.39,15.88,24.97,2.49 };
	vector < double > JPL_transfer_3 = { 14.16,24.94,2.87,8.10,9.00,23.13,23.09,23.09,22.83,24.98,24.98,24.93,24.94,9.10,13.44,24.99,24.94,24.99,24.98,24.96 };
	vector < double > JPL_transfer_4 = { 23.96,6.48,16.72,23.97,23.95,23.95,23.96,23.99,23.94,23.96 };
	vector < double > JPL_transfer_5 = { 0.45,3.17,24.93,10.34,12.53,7.11,13.44,24.94,24.94,24.98,22.19,24.99,22.01 };
	vector < double > JPL_transfer_6 = { 24.91,0.30,18.39,3.08,20.24,24.96,24.85,24.97,0.28 };
	vector < double > JPL_transfer_7 = { 15.69,0.50,9.83,24.94,24.90,24.48,20.87,24.91,0.66 };
	vector < double > JPL_transfer_8 = { 10.03,24.00,2.83,24.99,24.99,24.96,21.19,24.98 };
	vector < double > JPL_transfer_9 = { 22.69,4.24,24.47,24.46,24.47,24.44,24.46,24.46,24.46,18.54,9.22 };
	vector < double > JPL_transfer_10 = { 0.81,11.59,7.66,1.11,17.46,6.47,20.47,24.47,3.99 };

	//Delta v
	vector < double > JPL_dv_1 = { 161.8,139.2,65.8,208.2,115.2,300.1,564.9,78.3,105.0,233.3,453.5,340.4,300.8 };
	vector < double > JPL_dv_2 = { 659.0,301.1,252.1,143.8,146.8,68.6,40.6,84.2,105.3,448.5,148.0 };
	vector < double > JPL_dv_3 = { 219.1,80.8,105.2,55.2,140.2,85.5,95.0,237.6,205.9,149.9,245.2,71.6,197.3,160.4,132.2,240.0,161.2,364.3,230.4,232.5 };
	vector < double > JPL_dv_4 = { 86.1,103.1,62.6,222.9,709.1,553.9,219.9,233.9,739.0,232.6 };
	vector < double > JPL_dv_5 = { 129.6,45.2,172.9,52.6,160.7,280.8,221.1,163.5,98.2,115.7,164.8,674.8,291.1 };
	vector < double > JPL_dv_6 = { 156.0,198.0,305.8,71.2,194.4,920.5,314.1,353.0,272.8 };
	vector < double > JPL_dv_7 = { 400.6,173.6,211.3,374.4,109.6,171.2,145.1,194.3,233.0 };
	vector < double > JPL_dv_8 = { 287.9,111.9,112.2,144.5,540.0,260.1,198.8,82.7 };
	vector < double > JPL_dv_9 = { 83.3,148.1,495.9,464.9,405.2,285.9,254.8,62.3,156.6,36.5,174.9 };
	vector < double > JPL_dv_10 = { 189.4,112.9,110.0,121.3,117.9,280.1,300.4,120.6,70.2 };

	vector<vector<double>> JPL_dv;
	JPL_dv.push_back(JPL_dv_1);
	JPL_dv.push_back(JPL_dv_2);
	JPL_dv.push_back(JPL_dv_3);
	JPL_dv.push_back(JPL_dv_4);
	JPL_dv.push_back(JPL_dv_5);
	JPL_dv.push_back(JPL_dv_6);
	JPL_dv.push_back(JPL_dv_7);
	JPL_dv.push_back(JPL_dv_8);
	JPL_dv.push_back(JPL_dv_9);
	JPL_dv.push_back(JPL_dv_10);

	//Ts_Tf_eachmission
	vector<vector<double>> Ts(10), Tf(10);
	vector<vector<double>> JPL_rendezvous, JPL_transfer;
	JPL_rendezvous.push_back(JPL_rendezvous_1);
	JPL_rendezvous.push_back(JPL_rendezvous_2);
	JPL_rendezvous.push_back(JPL_rendezvous_3);
	JPL_rendezvous.push_back(JPL_rendezvous_4);
	JPL_rendezvous.push_back(JPL_rendezvous_5);
	JPL_rendezvous.push_back(JPL_rendezvous_6);
	JPL_rendezvous.push_back(JPL_rendezvous_7);
	JPL_rendezvous.push_back(JPL_rendezvous_8);
	JPL_rendezvous.push_back(JPL_rendezvous_9);
	JPL_rendezvous.push_back(JPL_rendezvous_10);
	JPL_transfer.push_back(JPL_transfer_1);
	JPL_transfer.push_back(JPL_transfer_2);
	JPL_transfer.push_back(JPL_transfer_3);
	JPL_transfer.push_back(JPL_transfer_4);
	JPL_transfer.push_back(JPL_transfer_5);
	JPL_transfer.push_back(JPL_transfer_6);
	JPL_transfer.push_back(JPL_transfer_7);
	JPL_transfer.push_back(JPL_transfer_8);
	JPL_transfer.push_back(JPL_transfer_9);
	JPL_transfer.push_back(JPL_transfer_10);
	for (int mission = 0; mission < 10; mission++)
	{
		Ts[mission].push_back(start_MJD_eachmission[mission] + JPL_rendezvous[mission][0]);
		Tf[mission].push_back(Ts[mission].back() + JPL_transfer[mission][0]);
		for (int i = 1; i < JPL_rendezvous[mission].size() - 1; i++)
		{
			Ts[mission].push_back(Tf[mission].back() + JPL_rendezvous[mission][i]);
			Tf[mission].push_back(Ts[mission].back() + JPL_transfer[mission][i]);
		}
	}


	//validation
	vector<double> T_temp, dv_temp;
	int mission_number = 0;
	T_temp.push_back(start_MJD_eachmission[mission_number]);
	for (int i = 0; i < Ts[mission_number].size(); i++)
	{
		T_temp.push_back(Ts[mission_number][i]);
		T_temp.push_back(Tf[mission_number][i]);
	}

	Dv_All(T_temp, JPL_sequence[mission_number], dv_temp);
	double total_dv = 0;
	for (int i = 0; i < dv_temp.size(); i++)
	{
		total_dv += dv_temp[i];
	}


	//Calculate cost function
	double total_socre = 0.0;
	for (int mission_id = 0; mission_id < JPL_dv.size(); mission_id++)
	{
		//Calculate mass
		double m0;
		double mp = mp_calc(JPL_dv[mission_id], m0);
		total_socre += 55.0 + 2.0e-6 * (m0 - 2000.0) * (m0 - 2000.0);
	}

	Opt_info_gtoc9 JPL_opt_info, Nan_opt_info;
	JPL_opt_info.debris_squence = JPL_sequence;
	//JPL_opt_info.debris_squence = Nan_sequence;
	Nan_opt_info.debris_squence = Nan_sequence;

	vector<double> T_opt_JPL;
	for (int mission = 0; mission < 10; mission++)
	{
		T_opt_JPL.push_back(start_MJD_eachmission[mission]);
		for (int i = 0; i < Ts[mission].size(); i++)
		{
			T_opt_JPL.push_back(Ts[mission][i]);
			T_opt_JPL.push_back(Tf[mission][i]);
		}
	}

	//double test = time_optimization_NLOPT( T_opt_JPL, JPL_opt_info.debris_squence);

	//test single mission (based on JPL solution)

	vector<vector<double>> T_single_mission_test,dv_missions_test;
	double start = clock();
	double socre = 0.0;
	//Nan_opt_info.debris_squence[0].pop_back();
	//Nan_opt_info.debris_squence[0].pop_back();
	for(int i =0; i< 0;i++)
		socre = DP_optimization_all_mission(JPL_opt_info.debris_squence, T_single_mission_test, dv_missions_test);
	//socre = DP_optimization_total(Nan_opt_info.debris_squence, T_single_mission_test, dv_missions_test);
	double end= clock();
	cout << " Score is " << socre << ";   time is " << (end - start) / CLOCKS_PER_SEC / 100.0 << endl;


	for (int i = 0; i < 0; i++)
		socre = DP_optimization_all_mission(Nan_opt_info.debris_squence, T_single_mission_test, dv_missions_test);
	end = clock();
	cout << " Score is " << socre << ";   time is " << (end - start) / CLOCKS_PER_SEC / 100.0 << endl;

	//vector<double> T_test_allmission;
	//for (int i = 0; i < T_single_mission_test.size(); i++) T_test_allmission.insert(T_test_allmission.end(), 
	//	T_single_mission_test[i].begin(), T_single_mission_test[i].end());
	//out_gtoc9_index_m0(T_test_allmission, &Nan_opt_info);
	//time_optimization_NLOPT(T_test_allmission, Nan_opt_info.debris_squence);
	//out_gtoc9_index_m0(T_test_allmission, &Nan_opt_info);

	vector < int > Test_1 = { 23, 113, 79, 22, 15, 50, 16, 118, 114, 38, 95 };
	vector < int > Test_2 = { 19, 48, 122, 12, 63, 7, 61, 107, 82, 45, 41, 39 };
	vector < int > Test_3 = { 72, 51, 69, 10, 73, 66, 28, 52, 64, 90, 58, 53, 56, 87, 111 ,35 };
	vector < int > Test_4 = { 108, 2 ,6, 0 ,18, 47 ,75 ,104, 26, 24 ,120 ,37 ,94 };
	vector < int > Test_5 = { 84, 59 ,11 ,1 ,98 ,40 ,62, 36 ,54 ,17 ,91 ,102 ,121, 76, 20 ,27, 80, 8 };
	vector < int > Test_6 = { 101, 78, 43, 74 ,5 ,103 };
	vector < int > Test_7 = { 49 ,46, 70, 9 ,105, 93, 109, 117, 55, 32, 119 };
	vector < int > Test_8 = { 86, 100, 34, 30, 65, 110, 115, 96, 81, 4, 92, 97, 88 };
	vector < int > Test_9 = { 33, 99, 68, 116, 71, 14, 106, 89, 60, 83 };
	vector < int > Test_10 = { 44, 67 ,42 ,57 ,3 ,21, 85 ,13 ,31, 77 ,25 ,29 ,112 };

	//vector < int > Test_1 = { 23, 50, 57, 96, 55, 121, 79, 113, 117, 25 };
	//vector < int > Test_2 = { 19, 45, 107, 61, 85, 41, 82, 115, 7, 43, 63, 21, 77, 88, 8 };
	//vector < int > Test_3 = { 72, 69, 51, 10, 73, 66, 28, 52, 64, 90, 29, 58, 87, 76 };
	//vector < int > Test_4 = { 108, 2, 6, 0, 122, 26, 75, 104, 37, 120, 24, 119, 32, 13 };
	//vector < int > Test_5 = { 84, 1, 98, 11, 59, 40, 47, 62, 36, 54, 17, 91, 102, 35, 112 };
	//vector < int > Test_6 = { 101, 60, 48, 103, 5, 39, 53, 74, 56, 38, 106 };
	//vector < int > Test_7 = { 49, 18, 9 ,46 ,105 ,93 ,70 ,31 ,109, 114, 4, 118 };
	//vector < int > Test_8 = { 86, 100, 34, 30, 65, 22, 92, 110, 81, 97 };
	//vector < int > Test_9 = { 33, 12, 68, 71, 89, 116, 14, 99, 20 ,27, 15, 95 };
	//vector < int > Test_10 = { 44, 42, 3, 67, 111, 78, 83, 80, 16, 94 };
	vector<vector<int>> Test_sequence;
	Test_sequence.push_back(Test_1);
	Test_sequence.push_back(Test_2);
	Test_sequence.push_back(Test_3);
	Test_sequence.push_back(Test_4);
	Test_sequence.push_back(Test_5);
	Test_sequence.push_back(Test_6);
	Test_sequence.push_back(Test_7);
	Test_sequence.push_back(Test_8);
	Test_sequence.push_back(Test_9);
	Test_sequence.push_back(Test_10);

	
	for(int i =0; i< Test_sequence.size();i++)
	{
		Test_sequence[i].pop_back();
	}
	Test_sequence[1].pop_back();
	Test_sequence[1].pop_back();

	opt_struct b;
	b.optimin = socre;
	b.dv = dv_missions_test;
	b.t = T_single_mission_test;
	//socre = Obj_gtoc9_aco(Test_sequence, &b);
	//socre = DP_optimization_total(Test_sequence, T_single_mission_test, dv_missions_test);
	cout << " Score is " << socre << ";   time is " << (end - start) / CLOCKS_PER_SEC / 100.0 << endl;

	auto test_sq1 = Test_sequence;
	double ab[2];

	//opt_struct a;
	//a.optimin = socre;
	//a.dv = dv_missions_test;
	//a.t = T_single_mission_test;
	//a.dv = b.dv;
	//a.t = b.t;
	//while(true)
	//{
	//	//double temp = localsearch_gtoc9_next(test_sq1, &a);
	//	//cout << " Score is " << temp << endl;
	//	if (temp < socre)
	//	{
	//		socre = temp;
	//		a.optimin = socre;
	//	}
	//	else
	//	{
	//		break;
	//	}
	//}
	

	cout << " Score is " << socre << ";   time is " << (end - start) / CLOCKS_PER_SEC / 100.0 << endl;
}


int main()
{
	load_input(); //读取空间碎片信息
	auto start_time = std::chrono::steady_clock::now();
	MultiTree multi_tree(128*20, 1, 20, 800, 10);
	multi_tree.Run();
	auto end_time = std::chrono::steady_clock::now();
	double duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / 1000.0;
	cout << "Time is " << duration << "s ; which is " << duration / 3600.0 << " hours" << endl;

	Test_GTOC9_results();
	return 1;
}

