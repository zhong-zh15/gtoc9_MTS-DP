#include "PSO.h"
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<omp.h>
#include<random>

#include "gtoc9_problem.h"


using namespace std;


/****************************************************************************
* 函数名   : PSO
* 功  能   : PSO优化函数
* 输  入   : X为优化变量，ConstX为参数， xbest为输入量输出为最终结果，fbest为输出指标，D为优化变量个数，Np为种群数，Itmax迭代总数
*            注意：所有变量必须归一化即x的范围是0-1
*            一般种群数为变量个数的10倍，迭代次数1000
****************************************************************************/
void PSO(double (*ObjFun)(double* X, const double* constX), const double* constX, double* xbest, double& fbest, int D, int Np, int ItMax, int ItOut,
	double OmegaMin, double OmegaMax, double C1Min, double C1Max, double C2Min, double C2Max, double Vmax)
{
	int ct;
	double Omega, C1, C2, val, GBVAL, rand0;
	double* PBVAL = new double[Np];
	double* GBPOS = new double[D];
	double** popvel = NULL;
	popvel = new double* [Np];
	for (int i = 0; i < Np; i++)
		popvel[i] = new double[D];

	double** pop = new double* [Np];
	for (int i = 0; i < Np; i++)
		pop[i] = new double[D];
	double** PBPOS = new double* [Np];
	for (int i = 0; i < Np; i++)
		PBPOS[i] = new double[D];

	for (int j = 0; j < Np; j++)
		PBVAL[j] = 1.0e10;
	GBVAL = 1.0e10;

	for (int i = 0; i < Np; i++)
	{
		for (int j = 0; j < D; j++)
		{
			rand0 =  realRand(0.0, 1.0);//第一个丢掉
			pop[i][j] =  realRand(0.0, 1.0);
			popvel[i][j] = 0.0;
		}
	}

	ct = 1;
	while (ct <= ItMax)
	{
		for (int i = 0; i < Np; i++)
		{
			val = ObjFun(pop[i], constX);
			if (val < PBVAL[i])
			{
				PBVAL[i] = val;
				for (int j = 0; j < D; j++)
					PBPOS[i][j] = pop[i][j];
			}
			{
				if (val < GBVAL)
				{
					GBVAL = val;
					for (int j = 0; j < D; j++)
						GBPOS[j] = pop[i][j];
				}
			}
		}
		Omega = OmegaMax - (OmegaMax - OmegaMin) / ItMax * ct;
		C1 = -(C1Max - C1Min) * ct / ItMax + C1Max;
		C2 = (C2Max - C2Min) * ct / ItMax + C2Min;
		for (int i = 0; i < Np; i++)
		{
			for (int j = 0; j < D; j++)
			{
				//rand0 =  realRand(0.0, 1.0);
				popvel[i][j] = Omega * popvel[i][j] + C1 * (PBPOS[i][j] - pop[i][j]) *  realRand(0.0, 1.0) + C2 * (GBPOS[j] - pop[i][j]) *  realRand(0.0, 1.0);
				if (popvel[i][j] > Vmax) popvel[i][j] = Vmax;//
				else if (popvel[i][j] < -Vmax) popvel[i][j] = -Vmax;//
				pop[i][j] += popvel[i][j];
				if ((pop[i][j] > 1.0) || (pop[i][j] < 0.0))
				{
					//rand0 =  realRand(0.0, 1.0);
					pop[i][j] =  realRand(0.0, 1.0);
					popvel[i][j] = 0.0;
				}
			}
		}
		for (int j = 0; j < D; j++)
			xbest[j] = GBPOS[j];
		fbest = GBVAL;
		if (ct % ItOut == 0)
		{
			cout << "No. of iteration=" << ct << endl;
			for (int i = 0; i < D; i++)
				cout << "xbest(" << i + 1 << ")=" << setprecision(15) << GBPOS[i] << endl;
			cout << "fbest=" << setprecision(15) << GBVAL << endl << endl;
		}
		ct = ct + 1;
	}

	delete[] PBVAL;
	delete[] GBPOS;
	for (int i = 0; i < Np; i++)
		delete[] popvel[i];
	delete[] popvel;
	for (int i = 0; i < Np; i++)
		delete[] pop[i];
	delete[] pop;
	for (int i = 0; i < Np; i++)
		delete[] PBPOS[i];
	delete[] PBPOS;
}


/****************************************************************************
* 函数名   : PSO
* 功  能   : PSO优化函数（重载版，输入为vector类型）
* 输  入   : X为优化变量，ConstX为参数， xbest为输入量输出为最终结果，fbest为输出指标，Np为种群数，Itmax迭代总数
*            注意：所有变量必须归一化即x的范围是0-1
*            一般种群数为变量个数的10倍，迭代次数1000
****************************************************************************/
void PSO(double(* ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* constX), void* constX,
	std::vector<double>& xbest, double& fbest, int Np, int ItMax, int ItOut, double OmegaMin, double OmegaMax, double C1Min,
	double C1Max, double C2Min, double C2Max, double Vmax)
{
	int D = xbest.size();
	int ct;
	double Omega, C1, C2, GBVAL, rand0;
	double* PBVAL = new double[Np];
	double* GBPOS = new double[D];
	double** popvel = NULL;
	popvel = new double* [Np];
	for (int i = 0; i < Np; i++)
		popvel[i] = new double[D];

	vector<vector<double>> pop(Np);
	for (int i = 0; i < pop.size(); i++) pop[i].resize(D);
	
	double** PBPOS = new double* [Np];
	for (int i = 0; i < Np; i++)
		PBPOS[i] = new double[D];

	for (int j = 0; j < Np; j++)
		PBVAL[j] = 1.0e10;
	GBVAL = 1.0e10;

	for (int i = 0; i < Np; i++)
	{
		for (int j = 0; j < D; j++)
		{
			rand0 = realRand(0.0, 1.0);//第一个丢掉
			//pop[i][j] = realRand(0.0, 1.0);
			popvel[i][j] = 0.0;
		}

		//test random
		Opt_info_gtoc9* JPL_opt_info = (Opt_info_gtoc9*)constX; //transfer to Opt_info_gtoc9       
		auto T_temp4opt = estimate_t_sequence_random(JPL_opt_info->debris_squence);
		inverse_map_input_gtoc9(T_temp4opt);
		for (int j = 0; j < D; j++) pop[i][j] = T_temp4opt[j];
	}

	ct = 1;
	while (ct <= ItMax)
	{
		
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < Np; i++)
		{
			vector<double> temp;
			double val = ObjFun(pop[i], temp, constX);
			if (val < PBVAL[i])
			{
				PBVAL[i] = val;
				for (int j = 0; j < D; j++)
				{
					PBPOS[i][j] = pop[i][j];
				}
			}

#pragma omp critical (PSO)
			{
				if (val < GBVAL)
				{
					GBVAL = val;
					for (int j = 0; j < D; j++)
					{
						GBPOS[j] = pop[i][j];
					}
				}
			}
		}

		Omega = OmegaMax - (OmegaMax - OmegaMin) / ItMax * ct;
		C1 = -(C1Max - C1Min) * ct / ItMax + C1Max;
		C2 = (C2Max - C2Min) * ct / ItMax + C2Min;
		for (int i = 0; i < Np; i++)
		{
			int counter = 0;
			for (int j = 0; j < D; j++)
			{
				popvel[i][j] = Omega * popvel[i][j] + C1 * (PBPOS[i][j] - pop[i][j]) * realRand(0.0, 1.0)
				                + C2 * (GBPOS[j] - pop[i][j]) * realRand(0.0, 1.0);
				if (popvel[i][j] > Vmax) popvel[i][j] = Vmax;//
				else if (popvel[i][j] < -Vmax) popvel[i][j] = -Vmax;//
				pop[i][j] += popvel[i][j];
				if ((pop[i][j] > 1.0) || (pop[i][j] < 0.0))
				{
					////rand0 =  realRand(0.0, 1.0);
					//pop[i][j] = realRand(0.0, 1.0);
					//popvel[i][j] = 0.0;
					counter = 1;
					break;
				}
			}

			if (counter == 1)
			{
				Opt_info_gtoc9* JPL_opt_info = (Opt_info_gtoc9*)constX; //transfer to Opt_info_gtoc9       
				auto T_temp4opt = estimate_t_sequence_random(JPL_opt_info->debris_squence);
				inverse_map_input_gtoc9(T_temp4opt);
				for (int j = 0; j < D; j++) pop[i][j] = T_temp4opt[j];
			}
		}
		for (int j = 0; j < D; j++)
			xbest[j] = GBPOS[j];
		fbest = GBVAL;
		if (ct % ItOut == 0)
		{
			cout << "No. of iteration=" << ct << endl;
			for (int i = 0; i < D; i++)
				cout << "xbest(" << i + 1 << ")=" << setprecision(15) << GBPOS[i] << endl;
			cout << "fbest=" << setprecision(15) << GBVAL << endl << endl;
		}
		ct = ct + 1;
	}

	for (int j = 0; j < D; j++)
		xbest[j] = GBPOS[j];
	
	fbest = GBVAL;

	
	delete[] PBVAL;
	delete[] GBPOS;
	for (int i = 0; i < Np; i++)
		delete[] popvel[i];
	delete[] popvel;

	for (int i = 0; i < Np; i++)
		delete[] PBPOS[i];
	delete[] PBPOS;
}

