/****************************************************************************
* Copyright (C), 2020-2031 Tsinghua University, School of Aerospace Engineering, LAD
* Author: Zhong Zhang
				zhong-zh19@mails.tsinghua.edu.cn
*               539977562@qq.com
* File: model.h
* Description: calculation J2 perturbed transfer cost
* (based on Simple ¦¤V Approximation for Optimization of Debris-to-Debris Transfers,
*  Hong-Xin Shen and Lorenzo Casalino,
*  Journal of Spacecraft and Rockets 2021 58:2, 575-580)
*
* Log:
*Version      Date        Author           Description
* 01        2022-03-15    Zhong Zhang       Create
****************************************************************************/

#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <math.h>
#include "Constant.h"


//global variable
extern double debris_data[123][8];   //derbris orbit data
extern double domega_debris[123];    
extern double dOmega_debris[123];

/****************************************************************************
* Function     : domega_init
* Description  : calculate omega drift rate of debris
*                input:
*					i: debris ID
*                ouput:
*					return value:  omega drift rate (unit: rad/s)
****************************************************************************/
double domega_init(int i);

/****************************************************************************
* Function     : dOmega_init
* Description  : calculate Omega drift rate of debris
*                input:
*					i: debris ID
*                ouput:
*					return value:  Omega drift rate (unit: rad/s)
****************************************************************************/
double dOmega_init(int i);

/****************************************************************************
* Function     : estimate_dv
* Description  : estimate delta v for debris-to-debris transfer, from debris_now to target starting from tnow;
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
double estimate_dv(int debris_now, int target, double tnow, double& Ts, double& Tf);
double estimate_dv(int debris_now, int target, double tnow, double& Ts, double& Tf, double end_epoch, int N);

/****************************************************************************
* Function     : Dv_ij
* Description  : calculate delta v for debris-to-debris transfer, from i to j, from tnow to tf;
*                (based on Simple ¦¤V Approximation for Optimization of Debris-to-Debris Transfers,
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
double Dv_ij(int i, int j, double ts, double tf);

/****************************************************************************
* Function     : Dv_All
* Description  : calculate delta v for all transfers
*                input:
*	                R: derbris sequence ID
*					T: time sequence (from 0, unit: Day)
*                ouput:
*				    dv_sequence: delta v sequence (unit: m/s)
*					return value: total delta v (unit: m/s)
****************************************************************************/
double Dv_All(const std::vector<double>& T, const std::vector<int>& R, std::vector<double>& dv_sequence);

#endif

