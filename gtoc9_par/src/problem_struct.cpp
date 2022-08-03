#include "problem_struct.h"

void Node_problem::operator=(const Node_problem& result)
{
	//ts_ = result.ts_;        //从上一节点出发的时刻
	//tf_ = result.tf_;        //到达下一节点的时刻
	next_debris_ = result.next_debris_;   //到达的空间碎片编号
	//dv_ = result.dv_;                     //此次全部转移的速度增量
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
