#include "problem_struct.h"

void Node_problem::operator=(const Node_problem& result)
{
	//ts_ = result.ts_;        //����һ�ڵ������ʱ��
	//tf_ = result.tf_;        //������һ�ڵ��ʱ��
	next_debris_ = result.next_debris_;   //����Ŀռ���Ƭ���
	//dv_ = result.dv_;                     //�˴�ȫ��ת�Ƶ��ٶ�����
}

void Solution_one::operator=(const Solution_one& a)
{
	debris_sequence_ = (a.debris_sequence_);
	t_sequence_ = a.t_sequence_;
	dv_sequence_ = a.dv_sequence_;
}

void Solution::operator=(const Solution& a)
{
	cost_ = a.cost_;
	total_removal_num_ = a.total_removal_num_;
	for (int i = 0; i < TreeNum; i++)
	{
		solution_[i] = a.solution_[i];
	}
}
