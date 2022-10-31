/****************************************************************************
* Copyright (C), 2020-2031 Tsinghua University, School of Aerospace Engineering, LAD
* Author: Zhong Zhang
				zhong-zh19@mails.tsinghua.edu.cn
*                 539977562@qq.com
* File: multitree_beam.cpp
* Description: The declaration class of MultiTree
*
* Log:
*Version      Date        Author           Description
* 01        2022-03-15    Zhong Zhang       Create
****************************************************************************/

#include "multitree_beam.h"

#include <cassert>
#include <numeric>

#include "DP.h"


/****************************************************************************
* Function     : G_Function
* Description  : calculate which tree in the TNC should be extended, only expand the tree with the lowest g score (not used in GTOC9)
*				 input: 
*					TNC
*				    counter
****************************************************************************/
void G_Function(const TNC& tnc, int& counter)
{
	double g_temp[TreeNum]; double g_min = 1.e20;  //Expand the tree with the smallest velocity increment
												  
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
* Function     : SortTNC
* Description  : sorting function: sort the h index from large to small according to the starting node number
*				 input:
*					TNC a
*				    TNC b
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
* Function     : EqualTNC
* Description  : determine whether two tnc sequences are equal
*                input:
*	                TNC a
*					TNC b
*                ouput:
*				    1 or 0
****************************************************************************/
inline bool MultiTree::EqualTNC(const TNC& a, const TNC& b)
{
	if (fabs(a.op_index_.cost_ - b.op_index_.cost_) < 1.0e-10)
	{
		//get two sequences
		Solution sa, sb;
		vector<int> a_sequence, b_sequence;
		//int afirst = 0, bfirst = 0;
		for (int i = 0; i < TreeNum; i++)
		{
			vector<Node*> solution_one_node_a, solution_one_node_b;
			a.tnc_[i]->return_node_sequence(solution_one_node_a);             //Returns a sequence of nodes, excluding the root node
			a.tnc_[i]->getback_problem(solution_one_node_a, sa.solution_[i]);//Backtracking function, returns the information of a node up to the root node
		
			a_sequence = sa.solution_[i].debris_sequence_;
			b.tnc_[i]->return_node_sequence(solution_one_node_b);             //Returns a sequence of nodes, excluding the root node
			b.tnc_[i]->getback_problem(solution_one_node_b, sb.solution_[i]);//Backtracking function, returns the information of a node up to the root node

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
* Function     : unique_remove
* Description  : deduplication function (input is sorted, output is the first W are deduplicated)
****************************************************************************/
void MultiTree::unique_remove(std::vector<TNC>& expandinglist)
{
	vector<TNC> new_expandinglist;
	new_expandinglist.reserve(W_);

	for (int i = 0; i<expandinglist.size()-1; i++) //The size will change after each deletion
	{
		if ( !EqualTNC(expandinglist[i],expandinglist[i+1])) //equal
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
* Function     : beam_discard
* Description  : sort function, and take the first W_
*                input:
*	                expandinglist
****************************************************************************/
void MultiTree::beam_discard(std::vector<TNC>& expandinglist)
{
	//按照h从大到小排序
	sort(expandinglist.begin(), expandinglist.end(), MultiTree::SortTNC);

	unique_remove(expandinglist);

}

/****************************************************************************
* Function     : init_max_dv
* Description  : the n with the smallest drift according to RAAN
*                input:
*	                init_expand_nodes: nodes at each level
*					visited: indicates whether the node is in the n columns. 0 is not, 1 is
*					n: column of indicator
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

	//Find the largest top n
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
* Function     : Initialize
* Description  : initialize the node to be expanded
*                input:
*	                expandinglist
****************************************************************************/
void MultiTree::Initialize(std::vector<TNC>& expandinglist)
{
	layer_ = TreeNum;
	vector<vector<Node*>> old_expand_nodes, new_expand_nodes;

	//the first level nodes of each tree are all expanded
	for (int i= 0; i< TreeNum; i++)
	{
		//generate all child nodes
		for (int j = 0; j < DebrisNum; j++) 
		{
			//Create new node
			Node_problem node_problem;
			node_problem.next_debris_ = j;
			Node* temp_node = new Node(multi_tree_[i].root_, node_problem);		
		}
		
	}


	////Select all in the first trees
	//old_expand_nodes.resize(DebrisNum);
	//for(int i = 0; i< DebrisNum; i++)
	//{
	//	old_expand_nodes[i].push_back(multi_tree_[0].root_->child_[i]);
	//}

	////The remaining trees are selected according to the maximum speed increment
	//for (int i = 1; i < TreeNum; i++)
	//{
	//	//For the current initial solution, get the next place to go
	//	for (int j =0; j< old_expand_nodes.size(); j++)
	//	{
	//		int visited[DebrisNum];
	//		init_max_dv(old_expand_nodes[j], visited, n_); //Pick the top n largest velocity increment transfers
	//		for (int k =0; k<DebrisNum; k++)
	//		{
	//			if(visited[k] == 1) //Indicates that this point needs to be expanded
	//			{
	//				vector<Node*> temp;
	//				for (auto temp_node : old_expand_nodes[j])
	//					temp.push_back(temp_node);            //include the previous node
	//				temp.push_back(multi_tree_[i].root_->child_[k]); //include this node
	//				new_expand_nodes.push_back(temp);
	//			}
	//		}
	//	}
	//	//Place the new extended point at the old point
	//	old_expand_nodes.clear();
	//	old_expand_nodes = new_expand_nodes;
	//	new_expand_nodes.clear();
	//}

	vector<Node*> temp;
	vector<vector<int>> s_vector = {
		{38 ,64 ,77 ,8 ,105 ,0 ,56 ,50 ,109 ,},
{38 ,64 ,45 ,43 ,112 ,36 ,48 ,50 ,109 ,},
{103 ,51 ,82 ,73 ,108 ,36 ,56 ,31 ,109 ,},
{103 ,51 ,77 ,73 ,112 ,89 ,48 ,97 ,12 ,},
{103 ,42 ,45 ,8 ,112 ,89 ,56 ,31 ,109 ,},
{38 ,42 ,77 ,73 ,108 ,36 ,56 ,50 ,116 ,},
{103 ,64 ,45 ,73 ,112 ,36 ,56 ,50 ,116 ,},
{118 ,51 ,45 ,73 ,105 ,89 ,48 ,31 ,12 ,},
{118 ,51 ,82 ,43 ,105 ,0 ,56 ,50 ,116 ,},
{103 ,42 ,82 ,43 ,108 ,89 ,48 ,50 ,12 ,},
{38 ,64 ,45 ,8 ,112 ,0 ,56 ,50 ,109 ,},
{38 ,51 ,45 ,73 ,112 ,36 ,56 ,31 ,109 ,},
{118 ,51 ,45 ,73 ,105 ,0 ,48 ,31 ,12 ,},
{118 ,64 ,45 ,8 ,112 ,36 ,48 ,50 ,109 ,},
{103 ,64 ,45 ,43 ,105 ,36 ,48 ,50 ,116 ,},
{38 ,64 ,77 ,73 ,112 ,36 ,48 ,50 ,12 ,},
{118 ,64 ,82 ,43 ,105 ,0 ,111 ,31 ,12 ,},
{103 ,51 ,45 ,8 ,112 ,0 ,111 ,50 ,116 ,},
{38 ,51 ,77 ,73 ,105 ,36 ,56 ,50 ,12 ,},
{118 ,51 ,77 ,8 ,105 ,0 ,56 ,50 ,109 ,},
{103 ,42 ,77 ,73 ,112 ,0 ,48 ,31 ,116 ,},
{118 ,42 ,82 ,73 ,105 ,89 ,48 ,97 ,12 ,},
{103 ,64 ,82 ,73 ,112 ,0 ,48 ,31 ,109 ,},
{118 ,51 ,82 ,73 ,108 ,0 ,48 ,97 ,109 ,},
{103 ,51 ,77 ,8 ,108 ,0 ,48 ,50 ,116 ,},
{103 ,64 ,82 ,73 ,105 ,0 ,56 ,97 ,109 ,},
{38 ,51 ,77 ,73 ,112 ,36 ,48 ,97 ,12 ,},
{118 ,64 ,45 ,43 ,105 ,0 ,56 ,31 ,12 ,},
{103 ,51 ,82 ,8 ,112 ,0 ,56 ,97 ,12 ,},
{103 ,51 ,77 ,8 ,105 ,36 ,48 ,31 ,116 ,},
{103 ,64 ,82 ,73 ,105 ,36 ,111 ,31 ,109 ,},
{38 ,64 ,82 ,43 ,112 ,0 ,56 ,31 ,109 ,},
{118 ,51 ,77 ,43 ,108 ,0 ,48 ,97 ,12 ,},
{103 ,111 ,7 ,8 ,75 ,31 ,101 ,97 ,102 ,},
{38 ,64 ,82 ,43 ,105 ,36 ,56 ,31 ,116 ,},
{103 ,42 ,82 ,43 ,105 ,89 ,48 ,31 ,109 ,},
{103 ,51 ,82 ,43 ,112 ,0 ,111 ,97 ,116 ,},
{38 ,64 ,82 ,43 ,112 ,36 ,56 ,31 ,109 ,},
{103 ,51 ,77 ,43 ,112 ,36 ,56 ,31 ,12 ,},
{38 ,111 ,77 ,73 ,75 ,31 ,101 ,91 ,59 ,},
{38 ,111 ,77 ,73 ,75 ,1 ,101 ,4 ,59 ,},
{38 ,42 ,77 ,8 ,82 ,36 ,33 ,97 ,102 ,},
{103 ,52 ,77 ,43 ,41 ,31 ,101 ,91 ,120 ,},
{38 ,42 ,7 ,8 ,75 ,31 ,33 ,91 ,102 ,},
{114 ,111 ,77 ,43 ,82 ,31 ,33 ,91 ,120 ,},
{38 ,52 ,7 ,43 ,41 ,36 ,78 ,97 ,120 ,},
{103 ,52 ,77 ,73 ,82 ,36 ,78 ,4 ,102 ,},
{103 ,52 ,77 ,43 ,41 ,1 ,33 ,91 ,120 ,},
{38 ,52 ,70 ,73 ,82 ,36 ,33 ,91 ,59 ,},
{103 ,42 ,7 ,43 ,41 ,1 ,78 ,4 ,59 ,},
	};




	for (int j = 0; j < s_vector.size(); j++)
	{
		for (int i = 0; i < TreeNum; i++)
		{
			auto s = s_vector[j];
			temp.push_back(multi_tree_[i].root_->child_[s[i]]);
		}
	}
	old_expand_nodes.clear();
	old_expand_nodes.push_back(temp);
	
	//Finally, complete the TNC
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

	//Clear the nodes that did not become TNC
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
* Function     : getback_problem
* Description  : backtracking function, returns the information of a node up to the root node
****************************************************************************/
inline void Node::getback_problem(vector<Node*>& solution_one_node, Solution_one & temp_solution_one)
{

	//The last one is the first node
	for (int i = solution_one_node.size() - 1; i > -1; i--)
	{
		temp_solution_one.debris_sequence_.push_back(solution_one_node[i]->problem_.next_debris_);

	}
}

/****************************************************************************
* Function     : return_node_sequence
* Description  : return clear sequence
****************************************************************************/
inline void Node::return_node_sequence(std::vector<Node*>& solution_one_node)
{
	for (Node* tempnode = this; tempnode->parent_; tempnode = tempnode->parent_)
	{
		solution_one_node.push_back(tempnode);
	}
}

/****************************************************************************
* Function     : Remove
* Description  : deletes a node and its parent, provided the node has no children; recursively
****************************************************************************/
inline void MultiTree::Remove(Node* a)
{
	Node* tempnode = a;

	while (true)
	{
		Node* tempfaternode = tempnode->parent_;

		if (tempnode->child_.size() == 0 && tempnode->inTNCcounter_ == 0 && tempnode->parent_) //The root node will not be deleted
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

/****************************************************************************
* Function     : Traverse
* Description  : traverse all nodes in the tree and return the node that needs to be deleted at the bottom
****************************************************************************/
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

/****************************************************************************
* Function     : delete_redundant_nodes
* Description  : delete all nodes in the tree that need to be deleted (redundant)
****************************************************************************/
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
* Function     : GetBack
* Description  : returns all information about all nodes in a TNC
****************************************************************************/
Solution TNC::GetBack()
{
	Solution temp;
	vector<vector<int>> debris_sequence(TreeNum);
	for (int i =0; i< TreeNum; i++)
	{
		std::vector<Node*> solution_one_node;
	    tnc_[i]->return_node_sequence(solution_one_node);             
		tnc_[i]->getback_problem(solution_one_node, temp.solution_[i]);//get debris sequence

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
* Function     : Calculate_op_index
* Description  : calculate index (each task optimizes the time separately, and the division time is equally divided)
****************************************************************************/
void TNC::Calculate_op_index()
{
	//TODO: just update with new results
	
	//optimal index
	int missionID = op_index_.expand_mission_ID;
	std::vector<Node*> solution_one_node;
	tnc_[missionID]->return_node_sequence(solution_one_node);               //get new sequence
	Solution_one temp;
	tnc_[missionID]->getback_problem(solution_one_node, temp);              //Backtracking to get the current result, only the sequence

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
* Function     : Calculate_op_index_change_time
* Description  : calculate index  (all tasks together optimize time)
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
		tnc_[i]->return_node_sequence(solution_one_node);               //get new sequence
		Solution_one temp;
		tnc_[i]->getback_problem(solution_one_node, temp);              //Backtracking to get the current result, only the sequence
		total_number += temp.debris_sequence_.size();                   //clean up

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
* Function     : ExpandNode
* Description  : expand a node
*                input:
*	                the node, and the observed sequence visited, the observed value is 1, otherwise it is 0
*					visited：the observed sequence, the observed is 1, otherwise it is 0
*                ouput:
*				    vector<Node*>: an extended list
****************************************************************************/
std::vector<Node*> Tree::ExpandNode(Node* node, const int* visited, const double* dv_nextdebris) 
//Expand the child nodes based on the problem, the nodes that have been repeated in this tree do not need to be repeated
{

	std::vector<Node*> expandnodes;                 //Expand child nodes

	std::vector<Node*> solution_one_node;
	node->return_node_sequence(solution_one_node);               //get the original node
	Solution_one node_debris_sequence;
	node->getback_problem(solution_one_node, node_debris_sequence);//Backtracking to get the current result, only the sequence

	int remove_num = 0;
	for (int i = 0; i < DebrisNum; i++)
	{
		if (visited[i] > 0) remove_num += 1;
	}

	for (int i = 0; i < DebrisNum; i++)
	{
		if (visited[i] > 0 || (remove_num <DebrisNum - TreeNum && dv_nextdebris[i] > dv_max_ ) ) continue; //already visited is greater than dv_max
		//if (visited[i] > 0 ) continue;
		
		auto ifchild = std::find_if(node->child_.begin(), node->child_.end(),
			[i](Node* a) {return a->key_ == i; }); //Determine if the point is within its child nodes

		if (ifchild == node->child_.end()) //child node does not exist
		{
			//Create new node
			Node_problem node_problem;
			node_problem.next_debris_ = i;
			//node_problem.dv_ = -1.0;
			Node* temp_node = new Node(node, node_problem);
			expandnodes.push_back(temp_node);           //put in new node
		}
		else //If it exists, put it in the existing child node
		{
			expandnodes.push_back(*ifchild);
		}
	}

	return expandnodes;
}

/****************************************************************************
* Function     : exclude_maxdv
* Description  : exclude nodes with too large velocity increments
****************************************************************************/
inline  void exclude_maxdv(Node* node, const int* visited, double* dv_debris, double tstart)
{
	std::vector<Node*> solution_one_node;
	node->return_node_sequence(solution_one_node);               //get the original node
	Solution_one node_debris_sequence;
	node->getback_problem(solution_one_node, node_debris_sequence);//Backtracking to get the current result, only the sequence
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
* Function     : Expansion_one_TNC
* Description  : extending a TNC
*                input:
*	                tnc: Former TNC
*                ouput:
*				    newTNCs: new extended TNCs
****************************************************************************/
void MultiTree::Expansion_one_TNC(const TNC& tnc, std::vector<TNC>& newTNCs)
{
	//The consideration is that all nodes can be expanded, but it is possible that the two tncs are the same after expansion
	newTNCs.clear();
	vector<TNC> delete_tnc;

	int visited[DebrisNum]{};                     //The observed sequence of this tnc, initialized to 0, and 1 after observation

	for (int j = 0; j < TreeNum; j++)             //Add 1 to all the visited nodes according to the target position
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

	for (int i = 0; i < TreeNum; i++) 
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
		exclude_maxdv(tnc.tnc_[i], visited, dv_debris,tnc.op_index_.t_mission[i][1]);     //Exclude fragments with large velocity increments

		omp_set_lock(&tnc.tnc_[i]->lock);
		vector<Node*> leafnode = multi_tree_[i].ExpandNode(tnc.tnc_[i], visited_temp, dv_debris); //child nodes extended by this node
		omp_unset_lock(&tnc.tnc_[i]->lock);
		
		for (int j = 0; j < leafnode.size(); j++) //Child node information extended by this node
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

	//After the expansion, get the indicator of the new TNC
	for (int i = 0; i < newTNCs.size(); i++)
	{
		if (layer_ >= DebrisNum - Change_time)
			newTNCs[i].Calculate_op_index_change_time();
		else
			newTNCs[i].Calculate_op_index();

	}

	//filter
	sort(newTNCs.begin(), newTNCs.end(), MultiTree::SortTNC);
	//Assign to expand list
	if (newTNCs.size() > b_)
	{
		newTNCs.erase(newTNCs.begin() + b_, newTNCs.end());
	}
}

/****************************************************************************
* Function     : Expansion
* Description  : get a new extension list based on extension information
****************************************************************************/
void MultiTree::Expansion(std::vector<TNC>& expandinglist)
{
	layer_++;                                             //current layer
	std::vector< std::vector<TNC>> new_childnode_private;               //local private variable
	
#pragma omp parallel
	{
		auto nthreads = omp_get_num_threads();
		auto id = omp_get_thread_num();

#pragma omp single
		{
			new_childnode_private.resize(nthreads);
		}

		//All nodes that need to be expanded
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
* Function     : Run
* Description  : running total function for multi-tree search
****************************************************************************/
void MultiTree::Run()
{
	for (int i = 0; i< TreeNum; i++)  multi_tree_[i].dv_max_ = dv_max_; 

	std::ofstream fout0("../output_result/result.txt");

	std::ofstream fout1("../output_result/All_results.txt");
	
	vector<TNC> expandinglist;
	Initialize( expandinglist);                         //Initialize the first expansion table
	
	while (expandinglist.size() > 0 && layer_ < DebrisNum) //非空
	{
		Expansion(expandinglist);

		Localsearch(expandinglist);

		beam_discard(expandinglist);                    //Sort and take the first W_

		delete_redundant_nodes();                       //delete redundant child nodes

		RecordBestResult(expandinglist,fout0);			//record the best

		if (layer_ == DebrisNum)
		{
			RecordAllResult(expandinglist, fout1);			         //record the best
		}
	}

	
}

/****************************************************************************
* Function     : RecordBestResult
* Description  : record the best results in each layer
****************************************************************************/
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


	
	//debris sequence
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
	//time sequence
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
	//speed increment sequence
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

/****************************************************************************
* Function     : RecordAllResult
* Description  : record all feasible results
****************************************************************************/
void MultiTree::RecordAllResult(std::vector<TNC>& expandinglist, std::ofstream& fout0)
{

   for(int i =0; i< expandinglist.size();i++)
   {

	   fout0 << "Removal total num and Total J:  " << expandinglist[i].op_index_.removal_num_<< " " << expandinglist[i].op_index_.cost_ << std::endl;

	   result_now_ = expandinglist[i].GetBack();


	   //std::cout << "Removal total num and Total J:  " << result_now_.total_removal_num_ << " " << result_now_.cost_ << std::endl;
	   fout0 << "Removal total num and Total J:  " << result_now_.total_removal_num_ << " " << result_now_.cost_ << std::endl;


   	
	   //debris sequence
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
	   //time sequence
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
	   //speed increment sequence
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

/****************************************************************************
* Function     : TNC_generation_by_nodes
* Description  : build TNC from sequence, time, speed increments
****************************************************************************/
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
				[next_debris](Node* a) {return a->key_ == next_debris; }); //Determine if the point is within its child nodes

			if (ifchild == node->child_.end()) //child node does not exist
			{
				//Create new node
				Node_problem node_problem;
				node_problem.next_debris_ = next_debris;
				Node* temp_node = new Node(node, node_problem);

				//omp_unset_lock(&node->lock);
				node = temp_node;
			}
			else //If it exists, put it in the existing child node
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

/****************************************************************************
* Function     : Localsearch
* Description  : local search of all TNCs in this layer
****************************************************************************/
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
			//if (remove_num >= DebrisNum - Last_mission) cout << "Local search end Score is " << a_all_end[i].optimin << endl;
		}
		n_all_end = X_all_end.size();
	}

	for (int i = 0; i < a_all_end.size(); i++)
	{
		//if (remove_num >= DebrisNum - Last_mission) cout << "Local search end Score is " << a_all_end[i].optimin << endl;
		expandinglist.push_back(TNC_generation_by_nodes(X_all_end[i], a_all_end[i]));
	}
}