/****************************************************************************
* Copyright (C), 2020-2031 Tsinghua University, School of Aerospace Engineering, LAD
* Author: Zhong Zhang
				zhong-zh19@mails.tsinghua.edu.cn
*                 539977562@qq.com
* File: OrbitFun.cpp
* Description: basic orbital function library for orbital parameters. This is developed by Fanghua Jiang
*
* Log:
*Version      Date        Author           Description
****************************************************************************/

#include "OrbitFun.h"

#include <math.h>       
#include <stdio.h>

#include "Constant.h"
#include "OrbitMath.h"

/****************************************************************************
* Function     : E2f
* Description  : calculate true anomaly based on eccentric anomaly and eccentricity
*                input:
*					E: eccentric anomaly (unit: rad): for hyperbolic orbits, refer to H, r=a(ecoshH-1)
*					e: eccentricity: circular orbit e=0, elliptical orbit 0<e<1,
*									 parabolic orbit e=1, hyperbolic orbit e>1
*                ouput:
*					return value: true anomaly (unit: rad)
*					flag: 1 or 0
****************************************************************************/
double E2f(int& flag, double E, double e)
{
	if (e < 0.0) { flag = 0; return E; }

	double f = 0.0;
	if (e >= 0.0 && e < 1.0)//circular and elliptical orbits
	{
		double E0 = fmod(E, D2PI);
		if (E0 > DPI)
			E0 -= D2PI;
		if (E0 < -DPI)
			E0 += D2PI;
		f = 2.0 * atan(sqrt((1.0 + e) / (1.0 - e)) * tan(0.5 * E0));
		f = f + E - E0;
	}
	else if (e > 1.0)//hyperbolic Orbit
		f = 2.0 * atan(sqrt((e + 1.0) / (e - 1.0)) * tanh(0.5 * E));
	else// (abs(e-1.0)<epsilon) Parabolic orbit
	{
		f = E;
	}
	flag = 1;
	return f;
}

/****************************************************************************
* Function     : E2M
* Description  : calculate mean anomaly based on eccentric anomaly and eccentricity
*                input:
*					E: eccentric anomaly (unit: rad), for hyperbolic orbits, refer to H
*					e: eccentricity: circular orbit e=0, elliptical orbit 0<e<1,
*									 parabolic orbit e=1, hyperbolic orbit e>1
*                ouput:
*					return value: mean anomaly (unit: rad), for hyperbolic orbits, refer to N
* 					flag: 1 or 0
****************************************************************************/
double E2M(int& flag, double E, double e)
{
	if (e < 0.0) { flag = 0; return E; }
	double M = 0.0;
	if (e >= 0.0 && e < 1.0)//circular and elliptical orbits
	{
		double E0 = fmod(E, D2PI);
		M = E0 - e * sin(E0);
		M = M + E - E0;
	}
	else if (e > 1.0)//hyperbolic Orbit
		M = e * sinh(E) - E;
	else//(abs(e-1.0)<epsilon) parabolic orbit
	{
		M = E;
	}
	flag = 1;
	return M;
}

/****************************************************************************
* Function     : f2E
* Description  : calculate eccentric anomaly based on ture anomaly and eccentricity
*                input:
*					f: true anomaly (unit: rad)
*					e: eccentricity: circular orbit e=0, elliptical orbit 0<e<1,
*									 parabolic orbit e=1, hyperbolic orbit e>1
*                ouput:
*					return value: eccentric anomaly, for hyperbolic orbits, refer to H, r=a(ecoshH-1)
* 					flag: 1 or 0
****************************************************************************/
double f2E(int& flag, double f, double e)
{
	if (e < 0.0) { flag = 0; return f; }

	double E = 0.0;
	if (e >= 0.0 && e < 1.0)//circular and elliptical orbits
	{
		double f0 = fmod(f, D2PI);
		if (f0 > DPI)
			f0 -= D2PI;
		if (f0 < -DPI)
			f0 += D2PI;
		E = 2.0 * atan(sqrt((1.0 - e) / (1.0 + e)) * tan(0.5 * f0));
		E += f - f0;
	}
	else if (e > 1.0)//hyperbolic Orbit
	{
		if (f > DPI - acos(1.0 / e) || f < -DPI + acos(1.0 / e))
		{
			flag = 0;
			return f;
		}
		else
			E = 2.0 * atanh(sqrt((e - 1.0) / (1.0 + e)) * tan(0.5 * f));
	}
	else//if(abs(e-1.0)<epsilon) parabolic orbit
	{
		E = f;
	}
	flag = 1;
	return E;
}

/****************************************************************************
* Function     : M2E
* Description  : calculate eccentric anomaly based on mean anomaly and eccentricity
*                input:
*					M: mean anomaly
*					e: eccentricity: circular orbit e=0, elliptical orbit 0<e<1,
*									 parabolic orbit e=1, hyperbolic orbit e>1
*					MaxIter: the maximum number of iterations, the default value is 100
*					epsilon: iterative absolute error, the default value is 1.0e-14
*                ouput:
*					return value: eccentric anomaly (unit: rad), for hyperbolic orbits, refer to H
* 					flag: 1 or 0
****************************************************************************/
double M2E(int& flag, double M, double e, int MaxIter, double epsilon)
{
	if (epsilon <= 0.0 || MaxIter < 1 || e < 0.0) { flag = 0; return M; }

	//For the iterative approach, see Chapter 2 of Solar System Dynamics, by Carl D. Murray and Stanley F. Dermott
	double E = 0.0, Minus = 0.0, DeMinus = 0.0, DeDeMinus = 0.0, DeDeDeMinus = 0.0, Delta1 = 0.0, Delta2 = 0.0, Delta3 = 0.0;
	int N = 0;
	if (e >= 0.0 && e < 1.0)//circular and elliptical orbits
	{
		double RM = fmod(M, D2PI);
		if (RM < 0.0)
			RM += D2PI;
		double sinRM = sin(RM);
		E = RM + 0.85 * e * Sign(sinRM);
		N = 0;
		Delta3 = 1.0;
		while (fabs(Delta3) >= epsilon && N < MaxIter)
		{
			Minus = E - e * sin(E) - RM;
			DeMinus = 1.0 - e * cos(E);
			DeDeMinus = e * sin(E);
			DeDeDeMinus = e * cos(E);
			Delta1 = -Minus / DeMinus;
			Delta2 = -Minus / (DeMinus + 0.5 * Delta1 * DeDeMinus);
			Delta3 = -Minus / (DeMinus + 0.5 * Delta2 * DeDeMinus + 1.0 / 6.0 * Delta2 * Delta2 * DeDeDeMinus);
			E = E + Delta3;
			N = N + 1;
		}
		E = E + M - RM;
	}
	else if (e > 1.0)//hyperbolic Orbit
	{
		E = asinh(M / e);
		Delta3 = 1.0;
		N = 0;
		while (fabs(Delta3) >= epsilon && N < MaxIter)
		{
			Minus = e * sinh(E) - E - M;
			DeMinus = e * cosh(E) - 1.0;
			DeDeMinus = e * sinh(E);
			DeDeDeMinus = e * cosh(E);
			Delta1 = -Minus / DeMinus;
			Delta2 = -Minus / (DeMinus + 0.5 * Delta1 * DeDeMinus);
			Delta3 = -Minus / (DeMinus + 0.5 * Delta2 * DeDeMinus + 1.0 / 6.0 * Delta2 * Delta2 * DeDeDeMinus);
			E = E + Delta3;
			N = N + 1;
		}
	}
	else //(abs(e-1.0)<epsilon) parabolic orbit
	{
		E = M;
	}
	if (((e >= 0.0 && e < 1.0) || (e > 1.0)) && fabs(Delta3) >= 5.0 * epsilon && N >= MaxIter)
	{
		flag = 0;
		return M;
	}
	flag = 1;
	return E;
}

/****************************************************************************
* Function     : f0dt2ft
* Description  : calculate the final true anomaly based on the initial true anomaly angle and evolution time
*                input:
*					f0: initial true anomaly (unit: rad)
*					dt: evolution time (unit: second)
*					a: semi-major axis (unit: meter), for parabolic orbits, a is the periapsis distance, i.e. p/2
*					e: eccentricity: circular orbit e=0, elliptical orbit 0<e<1,
*									 parabolic orbit e=1, hyperbolic orbit e>1
*					mu: the gravitational coefficient of the central planet 
*						(the product of the gravitational constant and the mass of the central planet), unit: m^3/s^2,
*					MaxIter: the maximum number of iterations, the default value is 100
*					epsilon: iterative absolute error, the default value is 1.0e-14
*                ouput:
*					return value: final true anomaly (unit: rad)
* 					flag: 1 or 0
****************************************************************************/
double f0dt2ft(int& flag, double f0, double dt, double a, double e, double mu, int MaxIter, double epsilon)
{
	if (mu <= 0.0 || MaxIter < 1 || a <= 0.0 || e < 0.0) { flag = 0; return f0; }

	double ft = 0.0;
	if ((e >= 0.0 && e < 1.0) || (e > 1.0))//circular, elliptical, hyperbolic orbit
	{
		double E = f2E(flag, f0, e);
		if (flag == 0) return f0;
		double M = E2M(flag, E, e);
		if (flag == 0) return f0;
		M += sqrt(mu / (a * a * a)) * dt;
		E = M2E(flag, M, e, MaxIter, epsilon);
		if (flag == 0) return f0;
		ft = E2f(flag, E, e);
		if (flag == 0) return f0;
	}
	else //(abs(e-1.0)<epsilon) parabolic orbit
	{
		if ((f0 < -DPI) || (f0 > DPI))
		{
			flag = 0; return f0;
		}
		else if (f0 > DPI || f0 < -DPI)
			ft = f0;
		else
		{
			double B = 0.75 * sqrt(2.0 * mu / (a * a * a)) * dt + 0.5 * tan(0.5 * f0) * ((tan(0.5 * f0)) * (tan(0.5 * f0)) + 3.0);
			double B1B = B + sqrt(1.0 + B * B);
			double tanv = 0.0;
			if (fabs(dt) < D2PI * sqrt((a * a * a) / mu) / 1000.0)//When the propulsion time is small
			{
				double A = pow(B1B, 2.0 / 3.0);
				tanv = 2.0 * A * B / (1.0 + (1.0 + A) * A);
			}
			else//When the propulsion time is not small
			{
				double temp = pow(B1B, 1.0 / 3.0);
				tanv = temp - 1.0 / temp;
			}
			ft = 2.0 * atan(tanv);
		}
	}
	flag = 1;
	return ft;
}

/****************************************************************************
* Function     : f0ft2dt
* Description  : calculate the evolution time based on the initial true anomaly and the final true anomaly
*                input:
*					f0: initial true anomaly (unit: rad)
*					ft: final true anomaly (unit: rad)
*					a: semi-major axis (unit: meter), for parabolic orbits, a is the periapsis distance, i.e. p/2
*					e: eccentricity: circular orbit e=0, elliptical orbit 0<e<1,
*									 parabolic orbit e=1, hyperbolic orbit e>1
*					mu: the gravitational coefficient of the central planet, the default value is Earth's gravitational coefficient 3.98600441500e+14,
*						(the product of the gravitational constant and the mass of the central planet), unit: m^3/s^2,
*                ouput:
*					return value: evolution time (unit: second)
* 					flag: 1 or 0
****************************************************************************/
double f0ft2dt(int& flag, double f0, double ft, double a, double e, double mu)
{
	flag = 0;
	if (mu <= 0.0 || a <= 0.0 || e < 0.0 ) //|| ft < f0)
		return 0.0;

	double dt = 0.0;
	if (e >= 1.0)
	{
		double maxangle = DPI - acos(1.0 / e);
		if ((Min(f0, ft) < -maxangle) || (Max(f0, ft) > maxangle))
		{		
			return 0.0;
		}
		else if (f0<-maxangle || f0>maxangle || ft<-maxangle || ft>maxangle)
		{
			dt = 1.0e308;
			return dt;
		}
	}

	double omega = sqrt(mu / (a * a * a));
	double delta = 0.0;
	if ((e >= 0.0 && e < 1.0) || (e > 1.0))
	{
		double E = f2E(flag, f0, e);
		double M0 = E2M(flag, E, e);
		E = f2E(flag, ft, e);
		double Mt = E2M(flag, E, e);
		if (flag == 0) return 0.0;
		delta = Mt - M0;
	}
	else// if(fabs(e-1.0)<epsilon)
	{
		double B1 = tan(0.5 * f0) * ((tan(0.5 * f0)) * (tan(0.5 * f0)) + 3.0);
		double B2 = tan(0.5 * ft) * ((tan(0.5 * ft)) * (tan(0.5 * ft)) + 3.0);
		delta = sqrt(2.0) / 3.0 * (B2 - B1);
	}
	dt = delta / omega;
	flag = 1;
	return dt;
}

/****************************************************************************
* Function     : coe2rv
* Description  : calculate the position and velocity components in the inertial Cartesian coordinate system
*                based on the classical orbital elements
*                input:
*					coe: classical orbital elements  coe[5]
*						 coe[0] semi-major axis a(unit: meter)
*						 coe[1] eccentricity e: circular orbit e=0, elliptical orbit 0<e<1,
*												parabolic orbit e=1, hyperbolic orbit e>1
*						 coe[2] orbital inclination i (unit: rad): 0<i<pi
*						 coe[3] right ascension of the ascending node Omega (unit: rad): when the orbital inclination is 0, it is meaningless, and its value can be set to 0
*						 coe[4] argument of periapsis omega (unit: rad): when the eccentricity is 0, it is meaningless, and its value can be set to 0
*						 coe[5] true anomaly f (unit: rad): when the eccentricity is 0, it is meaningless, and its value can be taken as omega+f
*					mu: the gravitational coefficient of the central planet
*                ouput:
*					rv: position and velocity components
* 					flag: 1 or 0
****************************************************************************/
void coe2rv(int& flag, double* rv, const double* coe, double mu = 3.98600441800e+14)
{
	flag = 0;

	double p = coe[0] * fabs(1.0 - coe[1] * coe[1]);//semi latus rectum: p

	double sini, cosi, sinO, cosO, sino, coso;
	sini = sin(coe[2]);
	cosi = cos(coe[2]);
	sinO = sin(coe[3]);
	cosO = cos(coe[3]);
	sino = sin(coe[4]);
	coso = cos(coe[4]);

	//Orbital plane normal unit vector, i.e. angular momentum unit vector
	double HVector[3] = { sini * sinO, -sini * cosO, cosi };

	//Eccentricity unit vector, also called Laplace vector
	double PVector[3] = { cosO * coso - sinO * sino * cosi, sinO * coso + cosO * sino * cosi, sino * sini };

	//semi latus rectum direction unit vector, PVector, QVector, HVector form a right-handed coordinate system
	// QVector=[-cosO*sino-sinO*coso*cosi;-sinO*sino+cosO*coso*cosi;coso*sini];
	double QVector[3];
	V_Cross(QVector, HVector, PVector);

	double r = 0.0;
	r = p / (1.0 + coe[1] * cos(coe[5]));

	for (int i = 0; i < 3; i++)
	{
		rv[i] = r * (cos(coe[5]) * PVector[i] + sin(coe[5]) * QVector[i]);
		rv[3 + i] = sqrt(mu / p) * (-sin(coe[5]) * PVector[i] + (cos(coe[5]) + coe[1]) * QVector[i]);
	}
	flag = 1;
	return;
}

/****************************************************************************
* Function     : rv2coe
* Description  : calculate the classical orbital elements based on the position and velocity components
*				 in the geocentric inertial Cartesian coordinate system
*                input:
*					RV: position and velocity components in the geocentric inertial Cartesian coordinate system
*					mu: the gravitational coefficient of the central planet
*                ouput:
*					coe: classical orbital elements
*						 coe[0] semi-major axis a (unit: meter)
*						 coe[1] eccentricity e: circular orbit e=0, elliptical orbit 0<e<1,
*												parabolic orbit e=1, hyperbolic orbit e>1
*						 coe[2] orbital inclination i (unit: rad): 0<i<pi
*				         coe[3] right ascension of the ascending node Omega (unit: rad): when the orbital inclination is 0, it is meaningless, and its value can be set to 0
*			             coe[4] argument of periapsis omega (unit: rad): when the eccentricity is 0, it is meaningless, and its value can be set to 0
*						 coe[5] true anomaly f (unit: rad): when the eccentricity is 0, it is meaningless, and its value can be taken as omega+f
* 					flag: 1 or 0
****************************************************************************/
void rv2coe(int& flag, double* coe, const double* RV, double mu = 3.98600441500e+14)
{
	int i;
	flag = 0;

	double R[3] = { RV[0], RV[1], RV[2] };
	double V[3] = { RV[3], RV[4], RV[5] };
	double radius = V_Norm2(R, 3);
	double velocity = V_Norm2(V, 3);

	double unitR[3];
	for (i = 0; i < 3; i++) unitR[i] = R[i] / radius;//radial unit vector
	double unitV[3];
	for (i = 0; i < 3; i++) unitV[i] = V[i] / velocity;//tangential unit vector
	double hvector[3];
	V_Cross(hvector, unitR, unitV);
	double h = radius * velocity * V_Norm2(hvector, 3);//angular momentum value
	
	double unith[3];
	for (i = 0; i < 3; i++) unith[i] = hvector[i] / V_Norm2(hvector, 3);//orbital plane normal unit vector
	//eccentricity vector
	double evector[3];
	V_Cross(evector, unitV, unith);
	for (i = 0; i < 3; i++) evector[i] = (velocity * h / mu) * evector[i] - unitR[i];
	coe[1] = V_Norm2(evector, 3);//eccentricity
	double p = h * h / mu;

		coe[0] = p / (fabs(1.0 - coe[1] * coe[1]));//semi-major axis a (unit: meter)
	bool judge = (coe[1] > 0.0);
	double unite[3] = { 0.0 };
	if (judge)
		for (i = 0; i < 3; i++) unite[i] = evector[i] / coe[1];//eccentricity unit vector
	coe[2] = acos(unith[2]);//orbital inclination i (unit: rad)

	double unitN[3] = { -unith[1], unith[0], 0.0 };

	double temp[3];

	if (V_Norm2(unitN, 3) == 0.0)
	{
		coe[3] = 0.0;
		if (!judge)
		{
			coe[4] = 0.0;       
			coe[5] = atan2(unitR[1] * unith[2], unitR[0]);
		}
		else
		{
			V_Cross(temp, unite, unitR);
			coe[4] = atan2(unite[1] * unith[2], unite[0]);        
			coe[5] = atan2(V_Dot(unith, temp, 3), V_Dot(unite, unitR, 3));
		}
	}
	else
	{
		V_Cross(temp, unitN, unitR);
		coe[3] = atan2(unith[0], -unith[1]);
		coe[5] = atan2(V_Dot(unith, temp, 3), V_Dot(unitN, unitR, 3));
		if (!judge)
		{
			coe[4] = 0.0;
		}
		else
		{
			V_Cross(temp, unitN, unite);
			coe[4] = atan2(V_Dot(unith, temp, 3), V_Dot(unite, unitN, 3));
			coe[5] = coe[5] - coe[4];
		}
	}
	//convert to [0,2pi)
	coe[3] = fmod(coe[3], D2PI);
	if (coe[3] < 0.0)
		coe[3] += D2PI;
	coe[4] = fmod(coe[4], D2PI);
	if (coe[4] < 0.0)
		coe[4] += D2PI;
	coe[5] = fmod(coe[5], D2PI);
	if (coe[1] >= 1.0)
	{
		if (coe[5] > DPI - acos(1.0 / coe[1]))
			coe[5] -= D2PI;
		else if (coe[5] < -DPI + acos(1.0 / coe[1]))
			coe[5] += D2PI;
	}
	flag = 1;
	return;
}

