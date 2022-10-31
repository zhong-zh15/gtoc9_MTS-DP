/****************************************************************************
* Copyright (C), 2020-2031 Tsinghua University, School of Aerospace Engineering, LAD
* Author: Zhong Zhang
				zhong-zh19@mails.tsinghua.edu.cn
*                 539977562@qq.com
* File: multitree_beam.h
* Description: The declaration class of MultiTree
*
* Log:
*Version      Date        Author           Description
* 01        2022-03-15    Zhong Zhang       Create
****************************************************************************/

#ifndef _BEAM_H
#define _BEAM_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iterator>
#include <omp.h>
#include <unordered_map>


#include "multitree_problem_struct.h"
#include "local_search.h"

/****************************************************************************
* Class        : Node
* Description  : used to store current satellite information
****************************************************************************/
class Node
{
public:
	Node* parent_;                       //parent node
	std::vector<Node*> child_;           //child node
	Node_problem problem_;               //the struct of the problem information
	unsigned int key_;                   //feature code, used to distinguish different child nodes under the parent node
	unsigned int inTNCcounter_;          //count the number of this node in the TNC

	omp_lock_t lock;
	//constructor
	Node()                               
	{
		parent_ = NULL; key_ = 0; inTNCcounter_ = 0; problem_.next_debris_ = -1; //problem_.dv_ = -1.0;
		omp_init_lock(&lock);
	}
	
	//constructor
	Node(Node* e, const Node_problem& a) 
	{    
		parent_ = e;                     //point to the parent node
		problem_ = a;

		e->child_.push_back(this);       //parent node points to child node
		
		this->key_ = a.next_debris_;     //assign a value to the feature code of this node
		inTNCcounter_ = 0;
		omp_init_lock(&lock);
	}

	//destructor: delete from bottom to top, that is, if the point has child nodes, it cannot be deleted
	~Node()                             
	{
		if (this->child_.size() == 0)
		{
			//find out where the node is (parent's child node) and delete it
			std::vector<Node*>& parent_child = (this->parent_)->child_; //temporary variables
			int value = this->key_;

			//omp_set_lock(&this->parent_->lock);
			
			auto iter = std::find_if(parent_child.begin(), parent_child.end(),
				[value](Node* temp) {return temp->key_ == value; });
			parent_child.erase(iter);
			
			//omp_unset_lock(&this->parent_->lock);
		}
		else
		{
			std::cout << "Child Node Exists But Is Deleted" << std::endl;
		}
		omp_destroy_lock(&lock);
	}

	void return_node_sequence(std::vector<Node*>& solution_one_node);       //returns a sequence of nodes, excluding the root node
	void getback_problem(vector<Node*>& solution_one_node, Solution_one & temp_solution_one);  //backtracking function, returns the information of a node up to the root node

};


/****************************************************************************
* Class        : Tree
* Description  : used to store current satellite information
****************************************************************************/
class Tree
{
public:
	Node* root_;
	double dv_max_;
	
	Tree()
	{
		dv_max_ = 0.0;
		root_ = new Node();
	}
	~Tree()
	{
		if (root_) delete root_;
	}

	std::vector<Node*> ExpandNode(Node*, const int* visited, const double* dv_nextdebris); //À©Õ¹º¯Êý
};

/****************************************************************************
* Class        : TNC
* Description  : contains multiple nodes, expressing the current state
****************************************************************************/
class TNC
{
public:
	Node *tnc_[TreeNum];               //array pointer, an array, which stores node pointers
	Optimization_index op_index_;

	TNC()
	{
		for (int i = 0; i< TreeNum; i++)
		tnc_[i] = nullptr;
	}
	
	TNC(Node* a[TreeNum])
	{
		for (int i = 0; i < TreeNum; i++)
			tnc_[i] = a[i];

		//the corresponding node TNCnumbe+1
		for (int i = 0; i< TreeNum; i++)
		{
			omp_set_lock(&a[i]->lock);
			a[i]->inTNCcounter_++;
			omp_unset_lock(&a[i]->lock);
		}
	}
	
	TNC(vector<Node*> a)
	{
		for (int i = 0; i < TreeNum; i++)
			tnc_[i] = a[i];

		//the corresponding node TNCnumbe+1
		for (int i = 0; i < TreeNum; i++)
		{
			omp_set_lock(&a[i]->lock);
			a[i]->inTNCcounter_++;
			omp_unset_lock(&a[i]->lock);
		}
	}
	
	//copy constructor
	TNC(const TNC & old)
	{
		op_index_ = old.op_index_;
		for (int i = 0; i< TreeNum; i++)
		{
			tnc_[i] = old.tnc_[i];
			//old.tnc_[i] = nullptr;
		}
		std::cout << "copy_construct " << std::endl;
	}

	//move constructor
	TNC(TNC && old) noexcept 
	{
		op_index_ = old.op_index_;
		for (int i = 0; i < TreeNum; i++)
		{
			tnc_[i] = old.tnc_[i];
			old.tnc_[i] = nullptr;
		}
	}

	TNC& operator=(const TNC& old) {
		
		op_index_ = old.op_index_;
		for (int i = 0; i < TreeNum; i++)
		{
			tnc_[i] = old.tnc_[i];
			//old.tnc_[i] = nullptr;
		}
		return *this;
	}

	TNC& operator=(TNC&& a)  noexcept
	{
		op_index_ = a.op_index_;
		for (int i = 0; i < TreeNum; i++)
		{
			tnc_[i] = a.tnc_[i];
			a.tnc_[i] = nullptr;
		}
		return *this;
	}
	
	~TNC()
	{

		if( tnc_[0])
		{
			//the corresponding node TNCnumbe-1
			for (int i = 0; i < TreeNum; i++)
			{
				omp_set_lock(&tnc_[i]->lock);
				tnc_[i]->inTNCcounter_--;
				omp_unset_lock(&tnc_[i]->lock);
			}
		}
	}

	Solution GetBack();  //get information on backtracking points for all nodes in the TNC

	void Calculate_op_index(); //calculate optimal metrics: number of cleans and total speed increment
	void Calculate_op_index_change_time();
};

/****************************************************************************
* Class        : MultiTree
* Description  : use of multiple tree frameworks
****************************************************************************/
class MultiTree
{
public:
	int layer_;                                                          //layer
	int W_;                                                              //beam width
	int n_;                                                              //initial node width
	int b_;                                                              //filter after a tnc extension
	double dv_max_;


	
	Tree multi_tree_[TreeNum];										     //store each tree

	Solution    result_now_;											 //the current generation optimal solution
	Solution    result_all_;											 //global optimal solution

	MultiTree(int beamwidth, int n, int b, double dvmax) :
		layer_(0), W_(beamwidth), n_(n), b_(b), dv_max_(dvmax){}   //constructor

	void Expansion_one_TNC(const TNC & tnc, std::vector<TNC> & newTNCs); //extending a TNC

	void Expansion(std::vector<TNC> & expandinglist);                    //extension function
	//void End_pso_optimize(vector<TNC>& expandinglist);

	static bool SortTNC(const TNC& a, const TNC& b);                     //sorted way
	static bool EqualTNC(const TNC& a, const TNC& b);
	void unique_remove(std::vector<TNC>& expandinglist);

	void beam_discard(std::vector<TNC>& expandinglist);                  //sort and take the first W_

	void Initialize(std::vector<TNC>& expandinglist);                    //initialize the first expansion table
	
	void Run();                                                          //run the main function

	void RecordBestResult(std::vector<TNC>& expandinglist,std::ofstream &fout0);                         //record the best solution

	void RecordAllResult(std::vector<TNC>& expandinglist, std::ofstream& fout0);
	
	TNC TNC_generation_by_nodes(vector<vector<int>>& X, opt_struct& para);


	void Localsearch(std::vector<TNC>& expandinglist);
	
	void Remove(Node* a);                                                //remove the node and its parent from the tree, provided there are no remaining child nodes
	
	void Traverse(Node* node, vector<Node*>& last_layer_nodes);          //traverse all nodes in the tree and return the node that needs to be deleted at the bottom

	void delete_redundant_nodes();                                       //delete all redundant nodes
};

//G_Function: calculate which tree in the TNC should be extended
void G_Function(const TNC& tnc, int& counter);

#endif




