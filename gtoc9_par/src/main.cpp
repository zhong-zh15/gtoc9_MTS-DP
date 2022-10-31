/****************************************************************************
* Copyright (C), 2020-2031 Tsinghua University, School of Aerospace Engineering, LAD
* Author: Zhong Zhang
				zhong-zh19@mails.tsinghua.edu.cn
*               539977562@qq.com
* File: main.cpp
* Description: the main program of the GTOC9 parallel computing example
*
* Log:
*Version      Date        Author           Description
* 01        2022-05-12    Zhong Zhang       Create
****************************************************************************/

#include "main.h"

#include <chrono>
#include "DP.h"
#include "gtoc9_problem.h"

//debris_data
//[0]: debris ID, [1]*: epoch, [2]: semi-major axis (meter), [3]: eccentricity, [4]: inclination (rad) 
//[5]: right ascension of the ascending node (rad), [6]: argument of perigee (rad), [7]: true anomaly (rad)
double debris_data[123][8];   //global debris information
const int debrisnum = 123;

double domega_debris[123];
double dOmega_debris[123];

/****************************************************************************
* Function     : load_input
* Description  : read space debris information and do preprocessing
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
 * 0 initial time
 * 1 a (unit：meter)
 * 2 e 
 * 3 i (unit：rad)
 * 4 W (unit：rad)
 * 5 omega (unit：rad)
 * 6 M (unit：rad)
 * 7 f (unit：rad)
 */

	for(int i = 0; i<123; i++)
	{
		double Re = req;// 6378137.0;	//Earth's equatorial radius (meter) 

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
		//debris_data
		//[0]: debris ID, [1]*: epoch, [2]: semi-major axis (meter), [3]: eccentricity, [4]: inclination (rad) 
        //[5]: right ascension of the ascending node (rad), [6]: argument of perigee (rad), [7]: true anomaly (rad)
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
		

	//the information needs to be converted into the initial moment, that is, MJD23467
	fin.close();

	
}


int main()
{
	load_input(); //read space debris information
	auto start_time = std::chrono::steady_clock::now();
	MultiTree multi_tree(128*20, 1, 20, 800);
	multi_tree.Run();
	auto end_time = std::chrono::steady_clock::now();
	double duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / 1000.0;
	cout << "Time is " << duration << "s ; which is " << duration / 3600.0 << " hours" << endl;

	return 1;
}

