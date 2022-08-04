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
* ��  ��   : Node
* ��  ��   : �ڵ��࣬���ڴ洢��ǰ������Ϣ
****************************************************************************/
class Node
{
public:
	Node* parent_;                       //���ڵ�
	std::vector<Node*> child_;           //�ӽڵ�
	Node_problem problem_;               //������Ϣ�Ľṹ��
	unsigned int key_;                   //�����룬�������ָ��ڵ��²�ͬ���ӽڵ�
	unsigned int inTNCcounter_;          //ͳ�Ƹýڵ���TNC�еĸ���

	omp_lock_t lock;
	//���캯��
	Node()                               
	{
		parent_ = NULL; key_ = 0; inTNCcounter_ = 0; problem_.next_debris_ = -1; //problem_.dv_ = -1.0;
		omp_init_lock(&lock);
	}
	
	//���캯��
	Node(Node* e, const Node_problem& a) 
	{    
		parent_ = e;                     //ָ�򸸽ڵ�
		problem_ = a;

		e->child_.push_back(this);       //��ָ���ӽڵ�
		
		this->key_ = a.next_debris_;     //�ýڵ�������븳ֵ
		inTNCcounter_ = 0;
		omp_init_lock(&lock);
	}

	//�������������¶���ɾ����������õ����ӽڵ㣬����ɾ��
	~Node()                             
	{
		if (this->child_.size() == 0)
		{
			//����ýڵ��� (���׵ĺ��ӽڵ㣩���ĸ�λ�ã�ɾ����
			std::vector<Node*>& parent_child = (this->parent_)->child_; //��ʱ����
			int value = this->key_;

			//omp_set_lock(&this->parent_->lock);
			
			auto iter = std::find_if(parent_child.begin(), parent_child.end(),
				[value](Node* temp) {return temp->key_ == value; });
			parent_child.erase(iter);
			
			//omp_unset_lock(&this->parent_->lock);
		}
		else
		{
			std::cout << "�����ӽڵ�ȴɾ��" << std::endl;
		}
		omp_destroy_lock(&lock);
	}

	void return_node_sequence(std::vector<Node*>& solution_one_node);            //���ؽڵ����У����������ڵ�
	void getback_problem(vector<Node*>& solution_one_node, Solution_one & temp_solution_one);        //���ݺ���������һ���ڵ�ֱ�����ڵ����Ϣ

};

/****************************************************************************
* ��  ��   : Tree
* ��  ��   : �����ϴ������нڵ�
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

	std::vector<Node*> ExpandNode(Node*, const int* visited, const double* dv_nextdebris); //��չ����
};

/****************************************************************************
* ��  ��   : TNC
* ��  ��   : ��������ڵ㣬��ﵱǰ״̬
****************************************************************************/
class TNC
{
public:
	Node *tnc_[TreeNum];               //����ָ�룬һ�����飬�����ŵ��ǽڵ�ָ��
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

		//��Ӧ�ڵ��TNCnumbe+1
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

		//��Ӧ�ڵ��TNCnumbe+1
		for (int i = 0; i < TreeNum; i++)
		{
			omp_set_lock(&a[i]->lock);
			a[i]->inTNCcounter_++;
			omp_unset_lock(&a[i]->lock);
		}
	}
	
	//�������캯��
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

	//�ƶ����캯��
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

		if( tnc_[0])  //��Ϊ��
		{
			//��Ӧ�ڵ��TNCnumbe-1
			for (int i = 0; i < TreeNum; i++)
			{
				omp_set_lock(&tnc_[i]->lock);
				tnc_[i]->inTNCcounter_--;
				omp_unset_lock(&tnc_[i]->lock);
			}
		}
	}

	

	Solution GetBack();  //��ȡTNC�����нڵ�Ļ��ݵ���Ϣ

	void Calculate_op_index(); //��������ָ��:������������ٶ�����
	void Calculate_op_index_change_time();
};

/****************************************************************************
* ��  ��   : MultiTree
* ��  ��   : ������ܵ�ʹ��
****************************************************************************/
class MultiTree
{
public:
	int layer_;                                                          //����
	int W_;                                                              //�������
	int n_;                                                              //��ʼ�ڵ���
	int b_;                                                              //һ��tnc��չ��ɸѡ
	double dv_max_;
	int range;

	
	Tree multi_tree_[TreeNum];										     //��Ÿ�����

	Solution    result_now_;											 //��ǰ�����Ž�
	Solution    result_all_;											 //ȫ�����Ž�

	MultiTree(int beamwidth, int n, int b, double dvmax, int range_ = 100) :
	layer_(0), W_(beamwidth), n_(n),b_(b),dv_max_(dvmax),range(range_){}   //���캯��

	void Expansion_one_TNC(const TNC & tnc, std::vector<TNC> & newTNCs); //��չһ��TNC

	void Expansion(std::vector<TNC> & expandinglist);                    //��չ����
	//void End_pso_optimize(vector<TNC>& expandinglist);

	static bool SortTNC(const TNC& a, const TNC& b);                     //��������
	static bool EqualTNC(const TNC& a, const TNC& b);
	void unique_remove(std::vector<TNC>& expandinglist);

	void beam_discard(std::vector<TNC>& expandinglist);                  //����ȡǰW_��

	void Initialize(std::vector<TNC>& expandinglist);                    //��ʼ���״���չ��
	
	void Run();                                                          //����������

	void RecordBestResult(std::vector<TNC>& expandinglist,std::ofstream &fout0);                         //��¼��ý�

	void RecordAllResult(std::vector<TNC>& expandinglist, std::ofstream& fout0);
	
	TNC TNC_generation_by_nodes(vector<vector<int>>& X, opt_struct& para);


	void Localsearch(std::vector<TNC>& expandinglist);
	
	void Remove(Node* a);                                                //������ɾ���ýڵ㼰�丸�ڵ㣬ǰ���ǲ����������ӽڵ�
	
	void Traverse(Node* node, vector<Node*>& last_layer_nodes);          //�����������нڵ㣬���صײ���Ҫɾ���Ľڵ�

	void delete_redundant_nodes();                                       //ɾ�����ж���ڵ�
};

//G_Function ����TNC���Ŀ���Ӧ�ñ���չ
void G_Function(const TNC& tnc, int& counter);

#endif




