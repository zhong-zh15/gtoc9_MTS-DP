/****************************************************************************
* Copyright (C), 2020-2031 Tsinghua University, School of Aerospace Engineering, LAD
* Author: Zhong Zhang
				zhong-zh19@mails.tsinghua.edu.cn
*                 539977562@qq.com
* File: problem_struct.h
* Description: classes on the structure of the solution
*
* Log:
*Version      Date        Author           Description
* 01        2021-05-12    Zhong Zhang       Create
****************************************************************************/

#include "multitree_problem_struct.h"

/****************************************************************************
* Function     : operator=
* Description  : equals overloading
****************************************************************************/
void Node_problem::operator=(const Node_problem& result)
{
	//ts_ = result.ts_;        //The moment from the previous node
	//tf_ = result.tf_;        //time to reach the next node
	next_debris_ = result.next_debris_;   //Arrived space debris ID
	//dv_ = result.dv_;                     //The speed increment for this entire transfer
}

/****************************************************************************
* Function     : operator=
* Description  : equals overloading
****************************************************************************/
void Solution_one::operator=(const Solution_one& a)
{
	debris_sequence_ = (a.debris_sequence_);
	t_sequence_ = a.t_sequence_;
	dv_sequence_ = a.dv_sequence_;
}

/****************************************************************************
* Function     : operator=
* Description  : equals overloading
****************************************************************************/
void Solution::operator=(const Solution& a)
{
	cost_ = a.cost_;
	total_removal_num_ = a.total_removal_num_;
	for (int i = 0; i < TreeNum; i++)
	{
		solution_[i] = a.solution_[i];
	}
}
