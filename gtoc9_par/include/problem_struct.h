/****************************************************************************
* Copyright (C), 2020-2031 �廪��ѧ���캽��ѧԺ����ѧ�����ʵ����
* ����: ���� zhong-zh19@mails.tsinghua.edu.cn
*            539977562@qq.com
* �ļ���: problem_struct.h
* ���ݼ��������ڽ�ṹ����
*
* �ļ���ʷ��
* �汾��     ����         ����       ˵��
* 01       2021-05-12    ����      �����ļ�
****************************************************************************/
#ifndef Problem_struct
#define Problem_struct

#include <vector>
#include <cstring>

#include "Constant.h"

/****************************************************************************
* �ṹ���� : Node_char
* ��  ��   : ��¼һ���ڵ�����������Ϣ
****************************************************************************/
class Node_problem
{
public:
	//double ts_;        //����һ�ڵ������ʱ��
	//double tf_;        //������һ�ڵ��ʱ��
	int next_debris_;    //����Ŀռ���Ƭ���
	//double dv_;          //�˴�ȫ��ת�Ƶ��ٶ�����

	Node_problem()
	{
		//ts_ = 0.0;
		//tf_ = 0.0;
		next_debris_ = -1;
		//dv_ = 0.0;
	}
	void operator =(const Node_problem& result);    //�Ⱥ�����
};



/****************************************************************************
* �ṹ���� : Optimization_index
* ��  ��   : �Ż�ָ��Ľṹ��
****************************************************************************/
class Optimization_index
{
public:
	int removal_num_;   //�ռ���Ƭ�������
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
* �ṹ���� : Solution_one
* ��  ��   : �������ǵĽ��
****************************************************************************/
class Solution_one
{
public:
	std::vector<int>    debris_sequence_;  //��������  N
	std::vector<double> t_sequence_;       //ʱ������  һ��2��N-1�����Ż�������������ʼʱ�̣��Լ�tstf
	std::vector<double> dv_sequence_;      //ÿ�������ٶ��������� N-1

	void operator =(const Solution_one& a); //��������
};


/****************************************************************************
* �ṹ���� : Solution
* ��  ��   : ���н��
****************************************************************************/
class Solution
{
public:
	Solution_one  solution_[TreeNum];      //���е����ǵĽ��
	double cost_;                      //�ܵ��ٶ�����
	int total_removal_num_;                //���������

	Solution():cost_(0.0), total_removal_num_(0){}
	
	Solution(const Solution& a) : cost_(a.cost_), total_removal_num_(a.total_removal_num_)
	{
		for (int i = 0; i< TreeNum; i++)
			solution_[i] = a.solution_[i];	
	}
	
	void operator =(const Solution& a);   //��������
};




#endif

