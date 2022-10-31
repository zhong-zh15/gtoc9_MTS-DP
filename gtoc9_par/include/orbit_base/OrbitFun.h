/****************************************************************************
* Copyright (C), 2020-2031 Tsinghua University, School of Aerospace Engineering, LAD
* Author: Zhong Zhang
				zhong-zh19@mails.tsinghua.edu.cn
*                 539977562@qq.com
* File: OrbitFun.h
* Description: basic orbital function library for orbital parameters. This is developed by Fanghua Jiang
*
* Log:
*Version      Date        Author           Description
****************************************************************************/

#ifndef _ORBITFUN_H_
#define _ORBITFUN_H_
#include <stdio.h>
#include"Constant.h"

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
double E2f(int& flag, double E, double e);

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
double E2M(int& flag, double E, double e);

/****************************************************************************
* Function     : f2E
* Description  : calculate eccentric anomaly based on ture anomaly and eccentricity
*                input:
*					f: true anomaly (unit: rad)
*					e: eccentricity: circular orbit e=0, elliptical orbit 0<e<1, 
*									 parabolic orbit e=1, hyperbolic orbit e>1
*                ouput:
*					return value: eccentric anomaly (unit: rad), for hyperbolic orbits, refer to H, r=a(ecoshH-1)
* 					flag: 1 or 0
****************************************************************************/
double f2E(int& flag, double f, double e);

/****************************************************************************
* Function     : M2E
* Description  : calculate eccentric anomaly based on mean anomaly and eccentricity
*                input:
*					M: mean anomaly (unit: rad)
*					e: eccentricity: circular orbit e=0, elliptical orbit 0<e<1, 
*									 parabolic orbit e=1, hyperbolic orbit e>1
*					MaxIter: the maximum number of iterations, the default value is 100
*					epsilon: iterative absolute error, the default value is 1.0e-14
*                ouput:
*					return value: eccentric anomaly (unit: rad), for hyperbolic orbits, refer to H
* 					flag: 1 or 0
****************************************************************************/
double M2E(int& flag, double M, double e, int MaxIter = 100, double epsilon = 1.0e-14);

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
double f0dt2ft(int& flag, double f0, double dt, double a, double e, double mu, int MaxIter = 100, double epsilon = 1.0e-14);

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
double f0ft2dt(int& flag, double f0, double ft, double a, double e, double mu = 3.98600441500e+14);

/****************************************************************************
* Function     : coe2rv
* Description  : calculate the position and velocity components in the inertial Cartesian coordinate system 
*                based on the classical orbital elements
*                input:
*					coe: classical orbital elements  coe[5]
*						 coe[0] semi-major axis a (unit: meter)
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
void coe2rv(int& flag, double* rv, const double* coe, double mu);

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
void rv2coe(int& flag, double* coe, const double* RV, double mu);

#endif