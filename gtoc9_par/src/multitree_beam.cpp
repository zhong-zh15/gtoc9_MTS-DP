#include "multitree_beam.h"

#include <cassert>
#include <numeric>

#include "DP.h"


// g function
// only expand the tree with the lowest g score (not used in GTOC9)
void G_Function(const TNC& tnc, int& counter)
{
	double g_temp[TreeNum]; double g_min = 1.e20;  //扩展速度增量最小的那颗树
												  
	for (int i = 0; i < TreeNum; i++)
	{
		g_temp[i] = tnc.op_index_.m_mission[i];  //test 1.mass() 2.time 3.dv
		
		std::vector<Node*> solution_one_node;
		tnc.tnc_[i]->return_node_sequence(solution_one_node);
		Solution_one temp;
		tnc.tnc_[i]->getback_problem(solution_one_node, temp);
		
		g_temp[i] = tnc.op_index_.dv_mission[i] / (temp.debris_sequence_.size());
		if (g_temp[i] < g_min) {
			g_min = g_temp[i]; counter = i;
		}
	}
}

/****************************************************************************
* 函数名   : SortNode(int a,int b)
* 功  能   : 排序函数：根据起始节点编号进行h指标从大到小排序
* 输 入    : 1.TNC a 2.TNC b
****************************************************************************/
inline bool MultiTree::SortTNC(const TNC& a, const TNC& b)
{
	double p = 1.0e-4;
	double h_a, h_b;
	double dt_a = 0.0, dt_b=0.0;

	if (a.op_index_.removal_num_ == 123)
	{
		h_a = a.op_index_.cost_;
		h_b = b.op_index_.cost_;
		return  h_a < h_b;
	}

	for (int i = 0; i < TreeNum; i++)
	{
		dt_a += a.op_index_.t_mission[i][1] - a.op_index_.t_mission[i][0];
		dt_b += b.op_index_.t_mission[i][1] - b.op_index_.t_mission[i][0];
	}
	h_a = a.op_index_.cost_ * pow( dt_a, p);
	h_b = b.op_index_.cost_ * pow( dt_b, p);

	return  h_a < h_b;
}

/****************************************************************************
* 函数名   : EqualTNC(const TNC& a, const TNC& b)
* 功  能   : 判断两个tnc序列是否相等
****************************************************************************/
inline bool MultiTree::EqualTNC(const TNC& a, const TNC& b)
{
	if (fabs(a.op_index_.cost_ - b.op_index_.cost_) < 1.0e-10)
	{
		//得到两个序列
		Solution sa, sb;
		vector<int> a_sequence, b_sequence;
		//int afirst = 0, bfirst = 0;
		for (int i = 0; i < TreeNum; i++)
		{
			vector<Node*> solution_one_node_a, solution_one_node_b;
			a.tnc_[i]->return_node_sequence(solution_one_node_a);             //返回节点序列，不包含根节点
			a.tnc_[i]->getback_problem(solution_one_node_a, sa.solution_[i]);//回溯函数，返回一个节点直至根节点的信息
		
			a_sequence = sa.solution_[i].debris_sequence_;
			b.tnc_[i]->return_node_sequence(solution_one_node_b);             //返回节点序列，不包含根节点
			b.tnc_[i]->getback_problem(solution_one_node_b, sb.solution_[i]);//回溯函数，返回一个节点直至根节点的信息
			
			b_sequence = sb.solution_[i].debris_sequence_;
			if (a_sequence != b_sequence)
			{
				return false;
			}
		}
		return true;
	}
	return  false;
}

/****************************************************************************
* 函数名   : void MultiTree::unique_remove(std::vector<TNC>& expandinglist)
* 功  能   : 去重函数(输入是已排序过的，输出是前W个是去重的)
****************************************************************************/
void MultiTree::unique_remove(std::vector<TNC>& expandinglist)
{
	vector<TNC> new_expandinglist;
	new_expandinglist.reserve(W_);

	for (int i = 0; i<expandinglist.size()-1; i++) //每次删除后大小都会变
	{
		if ( !EqualTNC(expandinglist[i],expandinglist[i+1])) //相等
		{
			new_expandinglist.push_back(std::move(expandinglist[i]));
		}
		if (new_expandinglist.size() >= W_) break;
	}

	expandinglist.clear();
	expandinglist.reserve(new_expandinglist.size());
	std::move(new_expandinglist.begin(), new_expandinglist.end(), back_inserter(expandinglist));
}
/****************************************************************************
* 函数名   : beam_discard(std::vector<TNC>& expandinglist)
* 功  能   : 排序函数，并取前W_个
****************************************************************************/
void MultiTree::beam_discard(std::vector<TNC>& expandinglist)
{
	//按照h从大到小排序
	sort(expandinglist.begin(), expandinglist.end(), MultiTree::SortTNC);

	unique_remove(expandinglist);

}
/****************************************************************************
* 函数名   : init_max_dv(vector<vector<Node*>> & init_expand_nodes, int tree_id, int *visited)
* 功  能   : 根据RAAN漂移最小的n个
* 输  入   : init_expand_nodes 每一层的节点，visited 表示该节点是否位于这n个之列 0为不在，1在, n为指标之列
****************************************************************************/
void init_max_dv(vector<Node*> & init_expand_node, int *visited, int n)
{
	for (int i = 0; i < DebrisNum; i++) visited[i] = 0;

	double visited_dv[DebrisNum];
	memset(visited_dv, 0.0, sizeof(double) * DebrisNum);
	
	for (int i = 0; i< init_expand_node.size(); i++)
	{
#pragma omp parallel for schedule(dynamic)
		for (int j = 0; j < DebrisNum; j++)
		{
			if (j == init_expand_node[i]->problem_.next_debris_)
			{
				visited_dv[j] = -1.e5;
			}
			else
			{
				double Ts, Tf;
				double dvtemp_min = 1.e10;
				for (double tnow = 0.0; tnow < 10.0; tnow += 2.5)
				{
					double dvtemp = estimate_dv(init_expand_node[i]->problem_.next_debris_, j, tnow, Ts, Tf);
					if (dvtemp < dvtemp_min)  dvtemp_min = dvtemp;
				}
				
				visited_dv[j] += dvtemp_min;
			}
		}
	}

	//找到最大的前n个
	double temp_dv[DebrisNum];
	memcpy(temp_dv, visited_dv, sizeof(double) * DebrisNum);
	std::sort(temp_dv, temp_dv + DebrisNum, std::greater<double>());
	double n_value = temp_dv[n-1];
	for (int i = 0; i< DebrisNum; i++)
	{
		if(visited_dv[i] > n_value - 1.0e-10)
		{
			visited[i] = 1;
		}
	}
	
}

/****************************************************************************
* 函数名   : Initialize(std::vector<TNC>& expandinglist)
* 功  能   : 初始化待扩展节点
****************************************************************************/
void MultiTree::Initialize(std::vector<TNC>& expandinglist)
{
	layer_ = TreeNum;
	//每棵树的第一层节点全部扩展
	for (int i= 0; i< TreeNum; i++)
	{
		//生成所有子节点
		for (int j = 0; j < DebrisNum; j++) 
		{
			//创建新节点
			Node_problem node_problem;
			node_problem.next_debris_ = j;
			Node* temp_node = new Node(multi_tree_[i].root_, node_problem);		
		}
		
	}

	vector<vector<Node*>> old_expand_nodes, new_expand_nodes;

	//第一颗树全部选取
	old_expand_nodes.resize(DebrisNum);
	for(int i = 0; i< DebrisNum; i++)
	{
		old_expand_nodes[i].push_back(multi_tree_[0].root_->child_[i]);
	}

	//其余树根据最大速度增量选
	for (int i = 1; i < TreeNum; i++)
	{
		//对当前初始解，获取下一个可去的地方
		for (int j =0; j< old_expand_nodes.size(); j++)
		{
			int visited[DebrisNum];
			init_max_dv(old_expand_nodes[j], visited, n_); //选取前n个最大的速度增量转移
			for (int k =0; k<DebrisNum; k++)
			{
				if(visited[k] == 1) //说明该点需要扩展
				{
					vector<Node*> temp;
					for (auto temp_node : old_expand_nodes[j])
						temp.push_back(temp_node);            //将前面的节点包含进来
					temp.push_back(multi_tree_[i].root_->child_[k]); //包含该点
					new_expand_nodes.push_back(temp);
				}
			}
		}
		//将新的扩展的点放置于旧点处
		old_expand_nodes.clear();
		old_expand_nodes = new_expand_nodes;
		new_expand_nodes.clear();
	}

	vector<Node*> temp;
	vector<int> s = { 103,111,7,8,75,31,101,97,102};
	for (int i = 0; i < TreeNum; i++)
	{
		temp.push_back(multi_tree_[i].root_->child_[s[i]]);
	}
	old_expand_nodes.clear();
	old_expand_nodes.push_back(temp);
	
	//最后将tnc完成
	for(int i=0; i < old_expand_nodes.size();i++)
	{
		TNC temp(old_expand_nodes[i]);
		for(int j =0; j< TreeNum;j++)
		{
			temp.op_index_.t_mission[j][0] = temp.op_index_.t_mission[j][1] = (2947.0+35.0)/TreeNum*j;
		}
		temp.op_index_.removal_num_ = TreeNum;
		expandinglist.emplace_back(std::move(temp));
	}

	//清除没有成为TNC的节点
	for (int i =1; i< TreeNum;i++)
	{
		while (true)
		{
			int j = 0;
			int sizechild = multi_tree_[i].root_->child_.size();
			for (; j < sizechild; j++)
			{
				if(multi_tree_[i].root_->child_[j]->inTNCcounter_ == 0 )
				{
					delete multi_tree_[i].root_->child_[j];
					break;
				}
			}
			if (j >= sizechild - 1)
				break;
		}
	}
	
}


/****************************************************************************
* 函数名   : Tree::getback(Node* node)
* 功  能   : 回溯函数，返回一个节点直至根节点的信息
****************************************************************************/
inline void Node::getback_problem(vector<Node*>& solution_one_node, Solution_one & temp_solution_one)
{

	//最后一个是最开始的节点
	for (int i = solution_one_node.size() - 1; i > -1; i--)
	{
		temp_solution_one.debris_sequence_.push_back(solution_one_node[i]->problem_.next_debris_);

	}
}

/****************************************************************************
* 函数名   : Node::return_node_sequence(std::vector<Node*>& solution_one_node)
* 功  能   : 返回清除序列
****************************************************************************/
inline void Node::return_node_sequence(std::vector<Node*>& solution_one_node)
{
	for (Node* tempnode = this; tempnode->parent_; tempnode = tempnode->parent_)
	{
		solution_one_node.push_back(tempnode);
	}
}


/****************************************************************************
* 函数名   : Remove()
* 功  能   : 删除一个节点及其父节点，前提是该节点没有子节点；采用递归方式
****************************************************************************/
inline void MultiTree::Remove(Node* a)
{
	Node* tempnode = a;

	while (true)
	{
		Node* tempfaternode = tempnode->parent_;

		if (tempnode->child_.size() == 0 && tempnode->inTNCcounter_ == 0 && tempnode->parent_) //根节点不会删除
		{
			delete tempnode;
		}
		else
		{
			break;
		}
		tempnode = tempfaternode;
	}
}

//遍历树中所有节点，返回底部需要删除的节点
void MultiTree::Traverse(Node* node, vector<Node*>& last_layer_nodes)
{
	if (node->child_.size() == 0 && node->inTNCcounter_ == 0)
	{
		last_layer_nodes.push_back(node);
		return;
	}

	for (int i = 0; i < node->child_.size(); i++)
	{
		Traverse(node->child_[i], last_layer_nodes);
	}
}

//删除树中所有需要删除（冗余）的节点
void MultiTree::delete_redundant_nodes()
{
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < TreeNum; i++)
	{
		vector<Node*> last_layer_nodes;
		Traverse(multi_tree_[i].root_, last_layer_nodes);
		for (int j = 0; j < last_layer_nodes.size(); j++)
		{
			Remove(last_layer_nodes[j]);
		}
	}
}

/****************************************************************************
* 函数名   : GetBack()
* 功  能   : 返回一个TNC中所有节点的全部信息
****************************************************************************/
Solution TNC::GetBack()
{
	Solution temp;
	vector<vector<int>> debris_sequence(TreeNum);
	for (int i =0; i< TreeNum; i++)
	{
		std::vector<Node*> solution_one_node;
	    tnc_[i]->return_node_sequence(solution_one_node);             
		tnc_[i]->getback_problem(solution_one_node, temp.solution_[i]);//得到碎片序列

		debris_sequence[i] = temp.solution_[i].debris_sequence_;
	}

	vector<vector<double>> time_sequence, dv_sequence;
	
	temp.cost_ = DP_optimization_all_mission(debris_sequence, time_sequence, dv_sequence, 1.0e10);

	if (time_sequence.size() > 0)
	{
		for (int i = 0; i < TreeNum; i++)
		{
			temp.solution_[i].debris_sequence_ = debris_sequence[i];
			temp.solution_[i].t_sequence_ = time_sequence[i];
			temp.solution_[i].dv_sequence_ = dv_sequence[i];
		}
	}

	temp.total_removal_num_ = op_index_.removal_num_;

	return temp;
}

/****************************************************************************
* 函数名   : Calculate_op_index()
* 功  能   : 计算指标（每个任务单独优化时间，划分时间是均分的）
****************************************************************************/
void TNC::Calculate_op_index()
{
	//TODO:只需更新新的结果
	
	//最优指标
	int missionID = op_index_.expand_mission_ID;
	std::vector<Node*> solution_one_node;
	tnc_[missionID]->return_node_sequence(solution_one_node);               //得到新序列
	Solution_one temp;
	tnc_[missionID]->getback_problem(solution_one_node, temp);              //回溯得到当前结果，只有序列

	double start_epoch, end_epoch;
	start_epoch = op_index_.t_mission[missionID][0];
	if (missionID == TreeNum - 1) end_epoch = 2947.0;
	else end_epoch = op_index_.t_mission[missionID+1][0] -35.0;

	double temp_end = (temp.debris_sequence_.size() - 1) * 23.5 + start_epoch;
	if (temp_end < end_epoch) end_epoch = temp_end;

	auto DP_result = DP_optimization_single_mission_min_m0(temp.debris_sequence_, start_epoch, end_epoch);

	vector<double> &time_sequence = DP_result.T_single_mission;
	vector<double>& dv_sequence = DP_result.dv_single_mission;
	if(time_sequence.size() > 0)
	{
		int i = missionID;
		{
			op_index_.t_mission[i][1] = time_sequence.back();
			op_index_.dv_mission[i] = accumulate(dv_sequence.begin(), dv_sequence.end(), 0.0);
			op_index_.m_mission[i] = DP_result.m0;
		}
	}

	op_index_.removal_num_ += 1;
	op_index_.cost_ = 0.0;
	for (int mission_id = 0; mission_id < TreeNum; mission_id++)
	{
		//Calculate mass
		double m0 = op_index_.m_mission[mission_id];
		op_index_.cost_ += 55.0 + 2.0e-6 * (m0 - 2000.0) * (m0 - 2000.0);
	}

}

/****************************************************************************
* 函数名   : Calculate_op_index_change_time()
* 功  能   : 计算指标（所有任务一起优化时间）
****************************************************************************/
void TNC::Calculate_op_index_change_time()
{
	//最优指标
	double total_dv = 0.0;
	int  total_number = 0;
	vector<vector<int>> debris_sequence(TreeNum);
	for (int i = 0; i < TreeNum; i++)
	{
		std::vector<Node*> solution_one_node;
		tnc_[i]->return_node_sequence(solution_one_node);               //得到新序列
		Solution_one temp;
		tnc_[i]->getback_problem(solution_one_node, temp);              //回溯得到当前结果，只有序列
		total_number += temp.debris_sequence_.size();                   //清理个数

		debris_sequence[i] = temp.debris_sequence_;
	}

	vector<vector<double>> time_sequence,dv_sequence;
	op_index_.cost_ = DP_optimization_all_mission(debris_sequence, time_sequence, dv_sequence,1.0e10);
	op_index_.removal_num_ = total_number;

	if (time_sequence.size() > 0)
	{
		for (int i = 0; i < TreeNum; i++)
		{
			op_index_.t_mission[i][0] = time_sequence[i][0];
			op_index_.t_mission[i][1] = time_sequence[i].back();
			op_index_.dv_mission[i] = accumulate(dv_sequence[i].begin(), dv_sequence[i].end(), 0.0);
			mp_calc(dv_sequence[i], op_index_.m_mission[i]);
		}
	}
}


/****************************************************************************
* 函数名   : vector<Node*> Tree::ExpandNode(Node* node, int* visited_debris)
* 功  能   : 扩展一个节点
* 输 入    : 该节点node，以及已经观测的序列visited，观测过的为1，否则为0
* 输    出 : vector<Node*> 一个扩展的列表
****************************************************************************/
std::vector<Node*> Tree::ExpandNode(Node* node, const int* visited, const double* dv_nextdebris) //基于问题扩展子节点,本棵树中已经重复的节点就不需要重复了
{

	std::vector<Node*> expandnodes;                 //扩展子节点

	std::vector<Node*> solution_one_node;
	node->return_node_sequence(solution_one_node);               //得到原有节点
	Solution_one node_debris_sequence;
	node->getback_problem(solution_one_node, node_debris_sequence);//回溯得到当前结果，只有序列

	int remove_num = 0;
	for (int i = 0; i < DebrisNum; i++)
	{
		if (visited[i] > 0) remove_num += 1;
	}

	for (int i = 0; i < DebrisNum; i++)
	{
		//if (visited[i] > 0 || dv_nextdebris[i] > dv_max_) continue; //已经访问过 大于dv_max

		if (visited[i] > 0 || (remove_num <DebrisNum - TreeNum && dv_nextdebris[i] > dv_max_ ) ) continue; //已经访问过 大于dv_max
		//if (visited[i] > 0 ) continue;
		
		auto ifchild = std::find_if(node->child_.begin(), node->child_.end(),
			[i](Node* a) {return a->key_ == i; }); //判断该点是否在其子节点内

		if (ifchild == node->child_.end()) //不存在子节点
		{
			//创建新节点
			Node_problem node_problem;
			node_problem.next_debris_ = i;
			//node_problem.dv_ = -1.0;
			Node* temp_node = new Node(node, node_problem);
			expandnodes.push_back(temp_node);           //放入新节点
		}
		else //存在即放入已有子节点
		{
			expandnodes.push_back(*ifchild);
		}
	}

	return expandnodes;
}


//排除掉速度增量太大的节点
inline  void exclude_maxdv(Node* node, const int* visited, double* dv_debris, double tstart)
{
	std::vector<Node*> solution_one_node;
	node->return_node_sequence(solution_one_node);               //得到原有节点
	Solution_one node_debris_sequence;
	node->getback_problem(solution_one_node, node_debris_sequence);//回溯得到当前结果，只有序列
	for (int j = 0; j < DebrisNum; j++)
	{
		if (visited[j] > 0)
		{
			dv_debris[j] = 1.0e10;
			continue;
		}
		double Ts, Tf;
		double dvtemp_min = 1.e10;
		//double tstart =  (node_debris_sequence.debris_sequence_.size() - 1) * 80.0;
		for (double tnow = tstart-150; tnow < tstart + 150; tnow += 10.0)
		{
			double dvtemp = estimate_dv(node_debris_sequence.debris_sequence_.back(), j, tnow, Ts, Tf);
			if (dvtemp < dvtemp_min)  dvtemp_min = dvtemp;
		}
		//if (node_debris_sequence.debris_sequence_.size() > 21 ) dvtemp_min = 0.0;
		
		dv_debris[j] = dvtemp_min;
	}

}



/****************************************************************************
* 函数名   : Expansion_one_TNC(const TNC& tnc, std::vector<TNC>& newTNCs)
* 功  能   : 扩展一个TNC
* 输 入    : const TNC& tnc 原TNC
* 输    出 : std::vector<TNC>& newTNCs  新扩展的tncs 
****************************************************************************/
void MultiTree::Expansion_one_TNC(const TNC& tnc, std::vector<TNC>& newTNCs)
{
	//考虑的是所有节点都可以扩展，但扩展后两个tnc是有可能相同的
	newTNCs.clear();
	vector<TNC> delete_tnc;

	int visited[DebrisNum]{};                     //该tnc的已观测序列，初始化为0，观测后为1

	for (int j = 0; j < TreeNum; j++)             //将访问过的所有节点按目标位置+1
	{
		std::vector<Node*> solution_one_node;
		tnc.tnc_[j]->return_node_sequence(solution_one_node);
		Solution_one temp;
		tnc.tnc_[j]->getback_problem(solution_one_node, temp);
		for (int k = 0; k < temp.debris_sequence_.size(); k++)
		{
			visited[temp.debris_sequence_[k]] ++;
		}
	}
	int counter = 0;


	G_Function(tnc, counter);

	for (int i = 0; i < TreeNum; i++) //按树扩展
	{
		//if (i != counter && layer_<DebrisNum-2*TreeNum) continue;

		int first_num = TreeNum;   //top first_num : if top 1 set 1
		first_num--;
		int visited_temp[DebrisNum]{};
		memcpy(visited_temp, visited, DebrisNum * sizeof(int));
		for (int derbis_id = 0; derbis_id < DebrisNum; derbis_id++)
		{
			if (visited[derbis_id] > 0) continue;
			
			vector<double> dv_min_temp(TreeNum);
			for (int tree_number = 0; tree_number < TreeNum; tree_number++)
			{
				int now_derbis = tnc.tnc_[tree_number]->key_;
				int target_debris = derbis_id;
				double now_time = tnc.op_index_.t_mission[tree_number][1];
				double dv_temp = Dv_ij(now_derbis, target_debris, now_time + 5.0, now_time + 30.0);
				dv_min_temp[tree_number] = dv_temp;
			}

			auto dv_min_temp_copy = dv_min_temp;
			sort(dv_min_temp_copy.begin(), dv_min_temp_copy.end());  //from low to high
			double dv_min_num = dv_min_temp_copy[first_num];

			if (dv_min_temp[i] > dv_min_num + 1.0e-10) // dv_min_tree expand, others not expanded
			{
				visited_temp[derbis_id] ++;
			}
		}

		
		double dv_debris[123];
		exclude_maxdv(tnc.tnc_[i], visited, dv_debris,tnc.op_index_.t_mission[i][1]);     //排除掉那些很大速度增量的碎片

		omp_set_lock(&tnc.tnc_[i]->lock);
		vector<Node*> leafnode = multi_tree_[i].ExpandNode(tnc.tnc_[i], visited_temp, dv_debris); //该节点扩展的子节点
		omp_unset_lock(&tnc.tnc_[i]->lock);
		
		for (int j = 0; j < leafnode.size(); j++) //该节点扩展的子节点信息
		{
			Node* temp[TreeNum];
			for (int k = 0; k < TreeNum; k++) temp[k] = tnc.tnc_[k];
			temp[i] = leafnode[j];
			TNC temp_tnc(temp);
			temp_tnc.op_index_ = tnc.op_index_;
			temp_tnc.op_index_.expand_mission_ID = i;
			newTNCs.push_back(std::move(temp_tnc));
		}
		
	}

	//扩展后，获得新TNC的指标
	for (int i = 0; i < newTNCs.size(); i++)
	{
		if (layer_ >= DebrisNum - Change_time)
			newTNCs[i].Calculate_op_index_change_time();
		else
			newTNCs[i].Calculate_op_index();

	}

	//先筛选
	sort(newTNCs.begin(), newTNCs.end(), MultiTree::SortTNC);
	//赋值到expand列表
	if (newTNCs.size() > b_)
	{
		newTNCs.erase(newTNCs.begin() + b_, newTNCs.end());
	}
}

/****************************************************************************
* 函数名   : Expansion(std::vector<TNC>& expandinglist) 扩展阶段函数
* 功  能   : 根据扩展信息得到新的扩展列表
****************************************************************************/
void MultiTree::Expansion(std::vector<TNC>& expandinglist)
{
	layer_++;                                             //当前层数
	std::vector< std::vector<TNC>> new_childnode_private;               //局部私有变量
	
#pragma omp parallel
	{
		auto nthreads = omp_get_num_threads();
		auto id = omp_get_thread_num();

#pragma omp single
		{
			new_childnode_private.resize(nthreads);
		}

		//所有需要扩展的节点
#pragma omp for schedule(dynamic) 
		for (int i = 0; i < expandinglist.size(); i++)
		{
			std::vector<TNC> newTNCs;
			Expansion_one_TNC(expandinglist[i], newTNCs);
			move(newTNCs.begin(), newTNCs.end(), back_inserter(new_childnode_private[id]));
		}
#pragma omp single
		{
			expandinglist.clear();
			for (auto& buffer : new_childnode_private) {
				move(buffer.begin(), buffer.end(), back_inserter(expandinglist));
				buffer.clear();
			}
		}
	}

}


/****************************************************************************
* 函数名   : Run()
* 功  能   : 多树搜索的运行总函数
****************************************************************************/
void MultiTree::Run()
{
	for (int i = 0; i< TreeNum; i++)  multi_tree_[i].dv_max_ = dv_max_; 

	std::ofstream fout0("../output_result/result.txt");

	std::ofstream fout1("../output_result/All_results.txt");
	
	vector<TNC> expandinglist;
	Initialize( expandinglist);                         //初始化首次扩展表
	
	while (expandinglist.size() > 0 && layer_ < DebrisNum) //非空
	{
		Expansion(expandinglist);

		Localsearch(expandinglist);

		beam_discard(expandinglist);                    //排序并取前W_个

		delete_redundant_nodes();                       //删除多余的子节点

		RecordBestResult(expandinglist,fout0);			//记录最好信息

		if (layer_ == DebrisNum)
		{
			RecordAllResult(expandinglist, fout1);			         //记录最好信息
		}
	}

	
}

//记录每一层中最好的结果
void MultiTree::RecordBestResult(std::vector<TNC>& expandinglist, std::ofstream& fout0)
{
	if (expandinglist.size() > 0)
	{
		std::cout << "Removal total num and Total dv:  " << expandinglist[0].op_index_.removal_num_ << " " << expandinglist[0].op_index_.cost_ << std::endl;
		fout0 << "Removal total num and Total dv:  " << expandinglist[0].op_index_.removal_num_ << " " << expandinglist[0].op_index_.cost_ << std::endl;

	}

	if (expandinglist.size()>0  ) 	result_now_ = expandinglist[0].GetBack();
	std::cout << "Removal total num and Total dv:  " << result_now_.total_removal_num_ << " " << result_now_.cost_ << std::endl;
	fout0 << "Removal total num and Total dv:  " << result_now_.total_removal_num_ << " " << result_now_.cost_ << std::endl;



	if (expandinglist.size() > 0)
	{
		vector<vector<int>> debris_sequence(TreeNum);
		Solution temp;
		int i = 0;
		for (int j = 0; j < TreeNum; j++)
		{
			std::vector<Node*> solution_one_node;
			expandinglist[i].tnc_[j]->return_node_sequence(solution_one_node);
			expandinglist[i].tnc_[j]->getback_problem(solution_one_node, temp.solution_[j]);//得到碎片序列
			debris_sequence[j] = temp.solution_[j].debris_sequence_;
		}
		vector<vector<double>> time_sequence, dv_sequence;
		double score_temp = DP_optimization_all_mission(debris_sequence, time_sequence, dv_sequence,1.0e10);
		std::cout << "Removal total num and Total J:  " << result_now_.total_removal_num_ << " " << score_temp << std::endl;
		fout0 << "Removal total num and Total J:  " << result_now_.total_removal_num_ << " " << score_temp << std::endl;
	}


	
	//碎片序列
	for (int i =0; i< TreeNum; i++) 
	{
		Solution_one &temp = result_now_.solution_[i];
		std::cout << "Removal num " << temp.debris_sequence_.size() << " : ";
		fout0 << "Removal num " << temp.debris_sequence_.size() << " : ";
		for (auto j: temp.debris_sequence_)
		{
			std::cout << j << ", ";
			fout0 << j << " ";
		}
		std::cout << std::endl;
		fout0 << std::endl;

	}
	//时间序列
	for (int i = 0; i < TreeNum; i++)
	{
		Solution_one& temp = result_now_.solution_[i];
		std::cout << "Epoch  " << temp.t_sequence_[0]<<"~" << temp.t_sequence_.back()<< " : ";
		fout0 << "Epoch  " << temp.t_sequence_[0] << "~" << temp.t_sequence_.back() << " : ";
		for (auto j : temp.t_sequence_)
		{
			std::cout << j << ", ";
			fout0 << j << " ";
		}
		std::cout << std::endl;
		fout0 << std::endl;
	}
	//速度增量序列
	for (int i = 0; i < TreeNum; i++)
	{
		Solution_one& temp = result_now_.solution_[i];
		double dv = 0.0;
		for (auto j : temp.dv_sequence_) dv += j;
		std::cout << "Dv " << dv << " : ";
		fout0 << "Dv " << dv << " : ";
		for (auto j : temp.dv_sequence_)
		{
			std::cout << j << ", ";
			fout0 << j << " ";
		}
		std::cout << std::endl;
		fout0 << std::endl;
	}
	
	std::cout << std::endl;
	fout0 << std::endl;

}

//记录所有可行的结果
void MultiTree::RecordAllResult(std::vector<TNC>& expandinglist, std::ofstream& fout0)
{

   for(int i =0; i< expandinglist.size();i++)
   {

	   fout0 << "Removal total num and Total J:  " << expandinglist[i].op_index_.removal_num_<< " " << expandinglist[i].op_index_.cost_ << std::endl;

	   result_now_ = expandinglist[i].GetBack();


	   //std::cout << "Removal total num and Total J:  " << result_now_.total_removal_num_ << " " << result_now_.cost_ << std::endl;
	   fout0 << "Removal total num and Total J:  " << result_now_.total_removal_num_ << " " << result_now_.cost_ << std::endl;


   	
	   //碎片序列
	   for (int i = 0; i < TreeNum; i++)
	   {
		   Solution_one& temp = result_now_.solution_[i];
		   //std::cout << "Removal num " << temp.debris_sequence_.size() << " : ";
		   fout0 << temp.debris_sequence_.size() << " ";
		   for (auto j : temp.debris_sequence_)
		   {
			   //std::cout << j << ", ";
			   fout0 << j << " ";
		   }
		   //std::cout << std::endl;
		   fout0 << std::endl;

	   }
	   //时间序列
	   for (int i = 0; i < TreeNum; i++)
	   {
		   Solution_one& temp = result_now_.solution_[i];
		   //std::cout << "Epoch  " << temp.t_sequence_[0] << "~" << temp.t_sequence_.back() << " : ";
		   //fout0 << "Epoch  " << temp.t_sequence_[0] << "~" << temp.t_sequence_.back() << " : ";

		   fout0.precision(10);
		   for (int j =0; j<  temp.t_sequence_.size(); j+=2)
		   {
			   double t_tmep = temp.t_sequence_[j];
			   //std::cout << j << ", ";
			   fout0 << t_tmep +23467.0 << " ";
		   }
		   //std::cout << std::endl;
		   fout0 << std::endl;
	   }
	   //速度增量序列
	   for (int i = 0; i < TreeNum; i++)
	   {
		   Solution_one& temp = result_now_.solution_[i];
		   double dv = 0.0;
		   for (auto j : temp.dv_sequence_) dv += j;
		   //std::cout << "Dv " << dv << " : ";
		   //fout0 << "Dv " << dv << " : ";
		   for (auto j : temp.dv_sequence_)
		   {
			   //std::cout << j << " ";
			   fout0 << j << " ";
		   }
		   //std::cout << std::endl;
		   fout0 << std::endl;
	   }

	   //std::cout << std::endl;
	   fout0 << std::endl;
   

   }

}


//根据序列、时刻、速度增量构建TNC
TNC  MultiTree::TNC_generation_by_nodes(vector<vector<int>> &X, opt_struct & para )
{
	Node* tnc_node[TreeNum];
	Optimization_index op_index_;
	op_index_.cost_ = para.optimin;
	int remove = 0;
	for (int mission = 0; mission < X.size(); mission++)
		remove += X[mission].size();
	op_index_.removal_num_ = remove;

	double t_mission[TreeNum];
	for (int j = 0; j < TreeNum; j++)
	{
		t_mission[j] = (2947.0 + 35.0) / TreeNum * j;
	}

	for(int mission =0; mission<X.size();mission++)
	{
		auto& debris_sequence = X[mission];
		auto& time_sequence = para.t[mission];
		auto& dv_sequence = para.dv[mission];

		op_index_.t_mission[mission][0] = t_mission[mission];
		op_index_.t_mission[mission][1] = time_sequence.back();
		op_index_.dv_mission[mission] = accumulate(dv_sequence.begin(), dv_sequence.end(), 0.0);
		mp_calc(dv_sequence, op_index_.m_mission[mission]);

		Node* node = multi_tree_[mission].root_;
		for(int i = 0; i<debris_sequence.size();i++)
		{
			int next_debris = debris_sequence[i];

			//omp_set_lock(&node->lock);

			auto ifchild = std::find_if(node->child_.begin(), node->child_.end(),
				[next_debris](Node* a) {return a->key_ == next_debris; }); //判断该点是否在其子节点内

			if (ifchild == node->child_.end()) //不存在子节点
			{
				//创建新节点
				Node_problem node_problem;
				node_problem.next_debris_ = next_debris;
				Node* temp_node = new Node(node, node_problem);

				//omp_unset_lock(&node->lock);
				node = temp_node;
			}
			else //存在即放入已有子节点
			{
				//omp_unset_lock(&node->lock);
				node = *ifchild;
			}
		}

		tnc_node[mission] = node;
	}

	if(op_index_.t_mission[0][0] > 1.0 )
	{
		cout << " Wrong!" << endl;
	}


	auto tnc = TNC(tnc_node);
	tnc.op_index_ = op_index_;
	return tnc;
}


//局部搜索该层的所有TNC
void MultiTree::Localsearch(std::vector<TNC>& expandinglist)
{
	int length = expandinglist.size();

	vector<vector<vector<int>>> X_all;
	vector<opt_struct> a_all;

	X_all.resize(length);
	a_all.resize(length);
	vector<vector<vector<int>>> X_all_end;
	vector<opt_struct> a_all_end;
	unordered_map<string, double> neighborhood;

	vector<vector<vector<int>>> X_all_temp;
	vector<opt_struct> a_all_temp;

	//get each tnc cost
#pragma omp parallel for schedule(dynamic)
	for (int num = 0; num < length; num++)
	{
		opt_struct& a = a_all[num];
		a.t.resize(TreeNum); a.dv.resize(TreeNum);
		a.optimin = expandinglist[num].op_index_.cost_;
		auto& tnc = expandinglist[num];
		vector<vector<int>> X(TreeNum);

		for (int i = 0; i < TreeNum; i++)
		{
			std::vector<Node*> solution_one_node;
			tnc.tnc_[i]->return_node_sequence(solution_one_node);
			Solution_one temp;
			tnc.tnc_[i]->getback_problem(solution_one_node, temp);
			X[i] = temp.debris_sequence_;

			if (!(tnc.op_index_.removal_num_ >= DebrisNum - Last_mission))
			{
				double start_epoch, end_epoch;
				start_epoch = tnc.op_index_.t_mission[i][0];
				if (i == TreeNum - 1) end_epoch = 2947.0;
				else end_epoch = tnc.op_index_.t_mission[i + 1][0] - 35.0;
				double temp_end = (temp.debris_sequence_.size() - 1) * 23.5 + start_epoch;
				if (temp_end < end_epoch) end_epoch = temp_end;
				auto DP_result = DP_optimization_single_mission_min_m0(X[i], start_epoch, end_epoch);
				a.t[i] = DP_result.T_single_mission;
				a.dv[i] = DP_result.dv_single_mission;
				double m0_temp; mp_calc(a.dv[i], m0_temp);
				a.opt_each_mission[i] = 55.0 + 2.0e-6 * (m0_temp - 2000.0) * (m0_temp - 2000.0);
			}
		}

		X_all[num] = X;
		if (tnc.op_index_.removal_num_ >= DebrisNum - Last_mission)
		{
			double score = DP_optimization_all_mission(X, a.t, a.dv, 10000.0);
#pragma omp critical
			{
				if (score < 10000.0 && a.t.size()>0 )
				{
					X_all_temp.push_back(X);
					a_all_temp.push_back(a);
				}
			}
		}
	}

	int remove_num = 0; for (int i = 0; i < X_all[0].size(); i++) remove_num += X_all[0][i].size();
	if (remove_num >= DebrisNum - Last_mission)
	{
		X_all.swap(X_all_temp); a_all.swap(a_all_temp); X_all_temp.clear(); a_all_temp.clear();
	}


	int n_all_end = 0;
	while (X_all.size() > 0)
	{
		double opt_max = -1.0e20;
		for (auto& a : a_all) if (a.optimin > opt_max) opt_max = a.optimin;
		opt_max = max(50.0, opt_max / TreeNum);
		opt_max = min(10000.0, opt_max);
		local_search_1layer(neighborhood, X_all, a_all, X_all_end, a_all_end, opt_max);
		for(int i= n_all_end; i< X_all_end.size();i++)
		{
			if (remove_num >= DebrisNum - Last_mission) cout << "Local search end Score is " << a_all_end[i].optimin << endl;
		}
		n_all_end = X_all_end.size();
	}

	for (int i = 0; i < a_all_end.size(); i++)
	{
		//if (remove_num >= DebrisNum - Last_mission) cout << "Local search end Score is " << a_all_end[i].optimin << endl;
		expandinglist.push_back(TNC_generation_by_nodes(X_all_end[i], a_all_end[i]));
	}
}