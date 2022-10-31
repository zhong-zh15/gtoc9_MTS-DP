/****************************************************************************
* Copyright (C), 2020-2031 Tsinghua University, School of Aerospace Engineering, LAD
* Author: Zhong Zhang
				zhong-zh19@mails.tsinghua.edu.cn
*               539977562@qq.com
* File: constant.h
* Description: constant table
*
* Log:
*Version      Date        Author           Description
* 01        2022-03-15    Zhong Zhang       Create
****************************************************************************/
#ifndef CONSTANT
#define CONSTANT
#include<math.h>

//global variable

const int DebrisNum = 123;
const int Mission_number = 9;
const int TreeNum = Mission_number;

const double Day2Second = 86400.0;      //seconds in a day (second)
const double mu = 398600.4418e9;	    //earth's gravitational coefficient (m^3/s^2)
const double req = 6378137.0;			//earth radius (meter)
const double J2 = 1.08262668e-3;		//J2 perturbation
const double MJD_Init= 23467.0;         //initial moment

//PI
const double pi = 3.1415926535897932384626433832795;	
const double DPI = 3.1415926535897932384626433832795;
const double D2PI = 6.283185307179586476925286766559;
const double D2R = 0.017453292519943295769236907684886;
const double R2D = 57.295779513082320876798154814105;

const int Last_mission = 0; //0 represents the last debris starts change time
const int Change_time = 0;

#endif

