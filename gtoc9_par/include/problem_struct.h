/****************************************************************************
* Copyright (C), 2020-2031 清华大学航天航空学院动力学与控制实验室
* 作者: 张众 zhong-zh19@mails.tsinghua.edu.cn
*            539977562@qq.com
* 文件名: problem_struct.h
* 内容简述：关于解结构的类
*
* 文件历史：
* 版本号     日期         作者       说明
* 01       2021-05-12    张众      创建文件
****************************************************************************/
#ifndef Problem_struct
#define Problem_struct

#include <vector>
#include <cstring>

#include "Constant.h"

/****************************************************************************
* 结构体名 : Node_char
* 功  能   : 记录一个节点关于问题的信息
****************************************************************************/
class Node_problem
{
public:
	//double ts_;        //从上一节点出发的时刻
	//double tf_;        //到达下一节点的时刻
	int next_debris_;    //到达的空间碎片编号
	//double dv_;          //此次全部转移的速度增量

	Node_problem()
	{
		//ts_ = 0.0;
		//tf_ = 0.0;
		next_debris_ = -1;
		//dv_ = 0.0;
	}
	void operator =(const Node_problem& result);    //等号重载
};



/****************************************************************************
* 结构体名 : Optimization_index
* 功  能   : 优化指标的结构体
****************************************************************************/
class Optimization_index
{
public:
	int removal_num_;   //空间碎片清理个数
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
* 结构体名 : Solution_one
* 功  能   : 单个卫星的结果
****************************************************************************/
class Solution_one
{
public:
	std::vector<int>    debris_sequence_;  //到访序列  N
	std::vector<double> t_sequence_;       //时间序列  一共2（N-1）个优化变量，包括初始时刻，以及tstf
	std::vector<double> dv_sequence_;      //每次所需速度增量序列 N-1

	void operator =(const Solution_one& a); //符号重载
};


/****************************************************************************
* 结构体名 : Solution
* 功  能   : 所有结果
****************************************************************************/
class Solution
{
public:
	Solution_one  solution_[TreeNum];      //所有的卫星的结果
	double cost_;                      //总的速度增量
	int total_removal_num_;                //总清理个数

	Solution():cost_(0.0), total_removal_num_(0){}
	
	Solution(const Solution& a) : cost_(a.cost_), total_removal_num_(a.total_removal_num_)
	{
		for (int i = 0; i< TreeNum; i++)
			solution_[i] = a.solution_[i];	
	}
	
	void operator =(const Solution& a);   //符号重载
};




#endif

