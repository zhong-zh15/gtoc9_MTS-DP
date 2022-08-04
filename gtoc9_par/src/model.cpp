/****************************************************************************
* Copyright (C), 2020-2031 Tsinghua University, School of Aerospace Engineering, LAD
* Author: Zhong Zhang
				zhong-zh19@mails.tsinghua.edu.cn
*               539977562@qq.com
* File: model.cpp
* Description: calcuation J2 perturbed transfer cost
* (based on Simple ΔV Approximation for Optimization of Debris-to-Debris Transfers, Hong-Xin Shen and Lorenzo Casalino, Journal of Spacecraft and Rockets 2021 58:2, 575-580)
*
* Log:
*Version      Date        Author           Description
* 01        2022-03-15    Zhong Zhang       Create
****************************************************************************/

#include "model.h"

/****************************************************************************
* Function     : domega_init
* Discription  : calculate omega drift rate of debris
*                input:
*					i: debris ID
*                ouput:
*					return value:  omega drift rate (unit: rad/s)
****************************************************************************/
double domega_init(int i)
{
	const double Re = req;// 6378137.0;													//m 地球赤道半径
//数据集debris_data
//0碎片编号 1* 2轨道半长轴(m) 3偏心率 4倾角(rad) 
//5升交点赤经(rad) 6近地点幅角(rad) 7真近点角(rad) 8雷达反射面积
	double a = debris_data[i][2];
	double e = debris_data[i][3];
	double inc = debris_data[i][4];
	double OMEGA = debris_data[i][5];
	double omega = debris_data[i][6];

	double p = a * (1 - e * e);
	double c2 = (Re / p) * (Re / p);
	double ci = cos(inc);
	double n = sqrt(mu / (a * a * a));
	double dOmega = -1.5 * J2 * c2 * n * ci;
	double dw = 0.75 * J2 * c2 * n * (5 * ci * ci - 1);
	domega_debris[i] = dw;
	return dw;
}

/****************************************************************************
* Function     : dOmega_init
* Discription  : calculate Omega drift rate of debris
*                input:
*					i: debris ID
*                ouput:
*					return value:  Omega drift rate (unit: rad/s)
****************************************************************************/
double dOmega_init(int i)
{
	const double Re = req;// 6378137.0;													//m 地球赤道半径
//数据集debris_data
//0碎片编号 1* 2轨道半长轴(m) 3偏心率 4倾角(rad) 
//5升交点赤经(rad) 6近地点幅角(rad) 7真近点角(rad) 8雷达反射面积
	double a = debris_data[i][2];
	double e = debris_data[i][3];
	double inc = debris_data[i][4];
	double OMEGA = debris_data[i][5];
	double omega = debris_data[i][6];

	double dOmega = -3.0 / 2.0 * J2 * (req / a) * (req / a) / (1.0 - e * e) / (1.0 - e * e) * sqrt(mu / a / a / a) * cos(inc);//碎片i的升交点赤经变化率
	dOmega_debris[i] = dOmega;
	return dOmega;
}

/****************************************************************************
* Function     : estimate_dv
* Discription  : estimate delta v for debris-to-debris transfer, from debris_now to target starting from tnow;
*                it will compute the suitable departure epoch (Ts) and arrival epoch (Tf)
*                input:
*					debris_now: debris now ID
*					target: target ID
*                   tnow: now time (from 0, unit: Day)
*                   end_epoch: the upper bound of the epoch (unit: Day,optional)
*					N: the number of left debris (optional)
*                ouput:
*					Ts: departure epoch (from 0, unit: Day)
*					Tf: arrival epoch (from 0, unit: Day)
*					return value: delta v (unit: m/s)
****************************************************************************/
double estimate_dv(int debris_now, int target, double tnow, double& Ts, double& Tf)
{
	double ts = (tnow  + 5.0) * Day2Second;						//转移起始时间ts设置为当前时间tnow的5天后
	double tf_low = ts + 1.0 * 86400.0;							//tf的下界，转移起始时间ts后的1天
	double tf_high = (tnow + 30.0)* Day2Second;					//tf的上界，当前时间tnow后的30天
	
	//计算转移至unvisited中碎片的速度增量估计值Estimatedv
	double a_s = debris_data[debris_now][2];
	double i_s = debris_data[debris_now][4];
	double Omega_s0 = debris_data[debris_now][5];
	double theta_s0 = debris_data[debris_now][6] + debris_data[debris_now][7];
	double a_f = debris_data[target][2];
	double i_f = debris_data[target][4];
	double Omega_f0 = debris_data[target][5];
	double theta_f0 = debris_data[target][6] + debris_data[target][7];
	
	double dOmega_s = -3.0 / 2.0 * J2 * (req / a_s) * (req / a_s) * sqrt(mu / a_s / a_s / a_s) * cos(i_s);
	double dOmega_f = -3.0 / 2.0 * J2 * (req / a_f) * (req / a_f) * sqrt(mu / a_f / a_f / a_f) * cos(i_f);

	//delta_Omega_low与delta_Omega_high表示升交点赤经差的变化范围
	double delta_Omega_low  = Omega_f0 + dOmega_f * tf_low  - Omega_s0 - dOmega_s * tf_low;		//tf_low 时刻，碎片debris_now与碎片unvisited[j]升交点赤经的差值
	double delta_Omega_high = Omega_f0 + dOmega_f * tf_high - Omega_s0 - dOmega_s * tf_high;	//tf_high时刻，碎片debris_now与碎片unvisited[j]升交点赤经的差值

	int delta_Omega_low_mod2pi = floor(delta_Omega_low / 2.0 / pi);
	int delta_Omega_high_mod2pi = floor(delta_Omega_high / 2.0 / pi);
	double tf_temp;
	if (delta_Omega_low_mod2pi != delta_Omega_high_mod2pi)
		//若delta_Omega_low与delta_Omega_high之间存在2pi的整数倍，则选择tf_temp使得升交点赤经差为2pi的整数倍
		tf_temp = (std::max(delta_Omega_low_mod2pi, delta_Omega_high_mod2pi) * 2.0 * pi
			- Omega_f0 + Omega_s0) / (dOmega_f - dOmega_s);
	else if (std::min(delta_Omega_low - delta_Omega_low_mod2pi * 2.0 * pi, delta_Omega_low_mod2pi * 2.0 * pi + 2.0 * pi - delta_Omega_low)
		< std::min(delta_Omega_high - delta_Omega_high_mod2pi * 2.0 * pi, delta_Omega_high_mod2pi * 2.0 * pi + 2.0 * pi - delta_Omega_high))
		//与delta_Omega_high相比，delta_Omega_low更靠近2pi的整数倍
		tf_temp = tf_low;
	else
		//与delta_Omega_low相比，delta_Omega_high更靠近2pi的整数倍
		tf_temp = tf_high;

	Ts = ts / Day2Second;
	Tf = tf_temp / Day2Second;
	return Dv_ij(debris_now, target, Ts, Tf);
}
double estimate_dv(int debris_now, int target, double tnow, double& Ts, double& Tf, double end_epoch, int N)
{
	double ts = (tnow + 5.0) * Day2Second;						//转移起始时间ts设置为当前时间tnow的5天后
	double tf_low = ts + 1.0 * 86400.0;							//tf的下界，转移起始时间ts后的1天
	double tf_high = (tnow + 30.0) * Day2Second;					//tf的上界，当前时间tnow后的30天

	
	double t_gap = (end_epoch - tnow) / N;
	if (t_gap < 30.0 ) tf_high = (tnow + t_gap) * Day2Second;					//tf的上界，当前时间tnow后的30天

	//计算转移至unvisited中碎片的速度增量估计值Estimatedv
	double a_s = debris_data[debris_now][2];
	double i_s = debris_data[debris_now][4];
	double Omega_s0 = debris_data[debris_now][5];
	double theta_s0 = debris_data[debris_now][6] + debris_data[debris_now][7];
	double a_f = debris_data[target][2];
	double i_f = debris_data[target][4];
	double Omega_f0 = debris_data[target][5];
	double theta_f0 = debris_data[target][6] + debris_data[target][7];

	double dOmega_s = -3.0 / 2.0 * J2 * (req / a_s) * (req / a_s) * sqrt(mu / a_s / a_s / a_s) * cos(i_s);
	double dOmega_f = -3.0 / 2.0 * J2 * (req / a_f) * (req / a_f) * sqrt(mu / a_f / a_f / a_f) * cos(i_f);

	//delta_Omega_low与delta_Omega_high表示升交点赤经差的变化范围
	double delta_Omega_low = Omega_f0 + dOmega_f * tf_low - Omega_s0 - dOmega_s * tf_low;		//tf_low 时刻，碎片debris_now与碎片unvisited[j]升交点赤经的差值
	double delta_Omega_high = Omega_f0 + dOmega_f * tf_high - Omega_s0 - dOmega_s * tf_high;	//tf_high时刻，碎片debris_now与碎片unvisited[j]升交点赤经的差值

	int delta_Omega_low_mod2pi = floor(delta_Omega_low / 2.0 / pi);
	int delta_Omega_high_mod2pi = floor(delta_Omega_high / 2.0 / pi);
	double tf_temp;
	if (delta_Omega_low_mod2pi != delta_Omega_high_mod2pi)
		//若delta_Omega_low与delta_Omega_high之间存在2pi的整数倍，则选择tf_temp使得升交点赤经差为2pi的整数倍
		tf_temp = (std::max(delta_Omega_low_mod2pi, delta_Omega_high_mod2pi) * 2.0 * pi
			- Omega_f0 + Omega_s0) / (dOmega_f - dOmega_s);
	//else if (std::min(delta_Omega_low - delta_Omega_low_mod2pi * 2.0 * pi, delta_Omega_low_mod2pi * 2.0 * pi + 2.0 * pi - delta_Omega_low)
	//	< std::min(delta_Omega_high - delta_Omega_high_mod2pi * 2.0 * pi, delta_Omega_high_mod2pi * 2.0 * pi + 2.0 * pi - delta_Omega_high))
	//	//与delta_Omega_high相比，delta_Omega_low更靠近2pi的整数倍
	//	tf_temp = tf_low;
	else
		//与delta_Omega_low相比，delta_Omega_high更靠近2pi的整数倍
		tf_temp = tf_high;

	Ts = ts / Day2Second;
	Tf = tf_temp / Day2Second;
	//return Dv_ij(debris_now, target, Ts, Tf);
	return 0.0;
}

/****************************************************************************
* Function     : Dv_ij
* Discription  : calculate delta v for debris-to-debris transfer, from i to j, from tnow to tf;
*                (based on Simple ΔV Approximation for Optimization of Debris-to-Debris Transfers,
*                 Hong-Xin Shen and Lorenzo Casalino,
*                 Journal of Spacecraft and Rockets 2021 58:2, 575-580)
*                input:
*					i: debris now ID
*					j: target ID
*					ts: departure epoch (from 0, unit: Day)
*					tf: arrival epoch (from 0, unit: Day)
*                ouput:
*					return value: delta v (unit: m/s)
****************************************************************************/
double Dv_ij(int i, int j, double ts, double tf)
{
	const double JD2S = 86400.0;													//天转秒
	const double muEarth = mu; // 3.98600441800e+14;								//m^3/s^2地球引力常数
	const double Re = req;// 6378137.0;												//m 地球赤道半径

	double mjd0 = ts * Day2Second;
	double mjdf = tf * Day2Second;
	//数据集debris_data
	//0碎片编号 1* 2轨道半长轴(m) 3偏心率 4倾角(rad) 
	//5升交点赤经(rad) 6近地点幅角(rad) 7真近点角(rad) 8雷达反射面积
	double a0 = debris_data[i][2];
	double af = debris_data[j][2];
	double e0 = debris_data[i][3];
	double ef = debris_data[j][3];
	double inc0 = debris_data[i][4];
	double incf = debris_data[j][4];
	//double dOmega_s = -3.0 / 2.0 * J2 * (req / a0) * (req / a0) * sqrt(mu / a0 / a0 / a0) * cos(inc0);//碎片i的升交点赤经变化率
	//double dOmega_f = -3.0 / 2.0 * J2 * (req / af) * (req / af) * sqrt(mu / af / af / af) * cos(incf);//碎片j的升交点赤经变化率
	// double Omega0 = debris_data[i][5] + dOmega_debris_notreal[i] * mjdf;  //tf时刻的值，0是出发碎片，f是到达碎片
	// double Omegaf = debris_data[j][5] + dOmega_debris_notreal[j] * mjdf;


	double Omega0 = debris_data[i][5] + dOmega_debris[i] * mjdf;  //tf时刻的值，0是出发碎片，f是到达碎片
	double Omegaf = debris_data[j][5] + dOmega_debris[j] * mjdf;
	double omega0 = debris_data[i][6] + domega_debris[i] * mjdf;
	double omegaf = debris_data[j][6] + domega_debris[j] * mjdf;
	double dt = (mjdf - mjd0);

	double a = (a0 + af) / 2;
	double inc = (inc0 + incf) / 2;
	double V = sqrt(muEarth / a);

	double delta_Omega = fmod(Omegaf - Omega0, 2 * pi);
	if (delta_Omega > pi)
		delta_Omega -= 2 * pi;
	else if (delta_Omega <= -pi)
		delta_Omega += 2 * pi;

	double sin_inc = sin(inc);
	double x = delta_Omega * sin_inc * V;
	double y = (af - a0) / 2. / a * V;
	double z = (incf - inc0) * V;

	//double Re_p0 = Re / a0 / (1 - e0 * e0);
	//double dOmega0 = -1.5 * J2 * Re_p0 * Re_p0 * sqrt(muEarth / a0 / a0 / a0) * cos(inc0);
	//double Re_pf = Re / af / (1 - ef * ef);
	//double dOmegaf = -1.5 * J2 * Re_pf * Re_pf * sqrt(muEarth / af / af / af) * cos(incf);
	double dOmega = (dOmega_debris[i] + dOmega_debris[j]) / 2.0;

	
	double m = 7 * dOmega * sin_inc * dt;
	//double n = dOmega * tan(inc) * sin(inc) * dt;
	double n = dOmega * sin_inc*sin_inc/cos(inc) * dt;

	double sx = (2. * x + m * y + n * z) / (4. + m * m + n * n);
	double sy = (2. * m * x - (4. + n * n) * y + m * n * z) / (8. + 2. * m * m + 2. * n * n);
	double sz = (2. * n * x + m * n * y - (4. + m * m) * z) / (8. + 2. * m * m + 2. * n * n);

	double delta_x = m * sy + n * sz;
	double dV1, dV2;

	// 偏心率修正
	double delta_ex = ef * cos(omegaf) - e0 * cos(omega0);
	double delta_ey = ef * sin(omegaf) - e0 * sin(omega0);
	//double delta_ex = ef * cos_omega[j] - e0 * cos_omega[i];
	//double delta_ey = ef * sin_omega[j] - e0 * sin_omega[i];
	double dVe = 0.5 * V * sqrt(delta_ex * delta_ex + delta_ey * delta_ey);
	dVe /= 2;
	dV1 = sqrt(sx * sx + sy * sy + sz * sz + dVe * dVe);
	dV2 = sqrt((x - sx - delta_x) * (x - sx - delta_x) + (y + sy) * (y + sy) + (z + sz) * (z + sz) + dVe * dVe);

	double dV;
	dV = dV1 + dV2;
	return dV;

}

/****************************************************************************
* Function     : Dv_All
* Discription  : calculate delta v for all transfers
*                input:
*	                R: derbris sequence ID
*					T: time sequence (from 0, unit: Day)
*                ouput:
*				    dv_sequence: delta v sequence (unit: m/s)
*					return value: total delta v (unit: m/s)
****************************************************************************/
double Dv_All(const std::vector<double>& T, const std::vector<int>& R, std::vector<double>& dv_sequence)
{

	int n = R.size();
	dv_sequence.resize(n - 1);
	for (int i = 0; i < n - 1; i++)
	{
		dv_sequence[i] = Dv_ij(R[i], R[i + 1], T[2 * i+1], T[2 * i + 2]);
	}
	return 0.0;
}