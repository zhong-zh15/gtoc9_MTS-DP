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

#ifndef Problem_struct
#define Problem_struct

#include <vector>
#include <cstring>

#include "Constant.h"

/****************************************************************************
* Class        : Node_char
* Description  : record a node's information about the problem
****************************************************************************/
class Node_problem
{
public:
	//double ts_;        //The initial moment from the previous node
	//double tf_;        //time to reach the next node
	int next_debris_;    //Arrived space debris ID
	//double dv_;          //The speed increment for this entire transfer

	Node_problem()
	{
		//ts_ = 0.0;
		//tf_ = 0.0;
		next_debris_ = -1;
		//dv_ = 0.0;
	}
	void operator =(const Node_problem& result);    //Equals overloading
};

/****************************************************************************
* Class        : Optimization_index
* Description  : class of optimization indicators
****************************************************************************/
class Optimization_index
{
public:
	int removal_num_;   //Number of space debris cleaned up
	double cost_;       //cost function
	int expand_mission_ID;
	double t_mission[TreeNum][2];
	double m_mission[TreeNum];
	double dv_mission[TreeNum];

	//int removal_num_each[TreeNum];
	int candidate_n_;

	Optimization_index()
	{
		removal_num_ = 0;
		cost_ = 0.0;
		expand_mission_ID = -1;
		candidate_n_ = 0;
		for(int i =0; i< TreeNum;i++)
		{
			t_mission[i][0] = 0.0;
			t_mission[i][1] = 0.0;
			m_mission[i] = 2000.0;
			dv_mission[i] = 0.0;
		}
	}
	Optimization_index(const Optimization_index & a): removal_num_(a.removal_num_),cost_(a.cost_), expand_mission_ID(a.expand_mission_ID),
		candidate_n_(a.candidate_n_)
	{

		memcpy(t_mission, a.t_mission, 2 * TreeNum * sizeof(double));
		memcpy(m_mission, a.m_mission, TreeNum * sizeof(double));
		memcpy(dv_mission, a.dv_mission, TreeNum * sizeof(double));
	}
	void operator =(const Optimization_index& a)
	{
		candidate_n_=(a.candidate_n_);
		expand_mission_ID=(a.expand_mission_ID);
		removal_num_ = (a.removal_num_);
		cost_ = (a.cost_);
	    memcpy(t_mission, a.t_mission, 2 * TreeNum * sizeof(double));
		memcpy(m_mission, a.m_mission, TreeNum * sizeof(double));
		memcpy(dv_mission, a.dv_mission, TreeNum * sizeof(double));
	}
};

/****************************************************************************
* Class        : Solution_one
* Description  : results for a single satellite
****************************************************************************/
class Solution_one
{
public:
	std::vector<int>    debris_sequence_;  //visiting sequence    N
	std::vector<double> t_sequence_;       //time sequence, a total of 2 (N-1) optimization variables, including the initial moment, and ts,tf
	std::vector<double> dv_sequence_;      //equence of each desired velocity increment    N-1

	void operator =(const Solution_one& a); //symbol overloading
};

/****************************************************************************
* Class        : Solution
* Description  : results for all solution
****************************************************************************/
class Solution
{
public:
	Solution_one  solution_[TreeNum];      //All satellite results
	double cost_;                      //total speed increment
	int total_removal_num_;                //total number of cleans

	Solution():cost_(0.0), total_removal_num_(0){}
	
	Solution(const Solution& a) : cost_(a.cost_), total_removal_num_(a.total_removal_num_)
	{
		for (int i = 0; i< TreeNum; i++)
			solution_[i] = a.solution_[i];	
	}
	
	void operator =(const Solution& a);   //symbol overloading
};

#endif

