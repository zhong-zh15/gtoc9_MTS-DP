#include "use_nlopt.h"


/****************************************************************************
* ������   : nloptmain()
* ��  ��   : ʹ��nlopt�Ż�ts��tf��ͬʱ��Ҫ����Լ����ʹ�÷�������
****************************************************************************/
void nloptmain(double (*ObjFun)(const std::vector<double>& x, std::vector<double>& grad, void* f_data), int* node_order,
	std::vector<double>& xbest, double& fbest, int ItMax)
{
	int num_variable = xbest.size();                                //ֻ�Ż�tf

	nlopt::opt opter(nlopt::LN_SBPLX, num_variable);				//�ֲ��Ż�  LN_SBPLX  LN_COBYLA LN_BOBYQA
		                                                                      //GN_DIRECT x  GD_STOGO x  GN_AGSx GN_CRS2_LM GN_ISRESx GN_ESCH
		                                                                      //Ч������ LN_SBPLX�� GN_CRS2_LM  �ȽϺ�

	opter.set_min_objective(ObjFun, node_order);					//ָ��
	
	std::vector<double> lb(num_variable), rb(num_variable), dx(num_variable);

	//opter.add_inequality_mconstraint(ObjFun, node_order, dx);
	
	fbest = 1.0e20;
	for (int i = 0; i < num_variable; i++)
	{
		lb[i] = 0.0;												//�½�
		rb[i] = 2952.0;												//�Ͻ�
		dx[i] = 1.0/86400.0*200.0;									//��ʼ����
	}
	//opter.set_initial_step(dx);										//���ó�ʼ����
	opter.set_lower_bounds(lb);										//�����½�
	opter.set_upper_bounds(rb);										//�����Ͻ�
	double tol = 1e-8;
	//opter.set_ftol_abs(tol);
	//opter.set_force_stop(tol);

	opter.set_maxeval(ItMax);										//�Ż��ô�����ֹͣ

	nlopt::result res = opter.optimize(xbest, fbest);				//�����Ż�

}


/****************************************************************************
* ������   : nlopt_tf()
* ��  ��   : ʹ��nlopt�Ż�tf(n-1���Ż�������tsֱ������)��ͬʱ��Ҫ����Լ����ʹ�÷��������Ͻ���½�
****************************************************************************/
double nlopt_tstf(const std::vector<double>& X, std::vector<double>& grad, void* f_data)
{
	int* node_order = (int*)f_data;
	int derbrisnum = (X.size() + 1 )/2 ;

	double penalty = 0.0;
	double tnow =  X[0]; //��ʼ����ʱ��
	double t_low = 5.0; //�ٴγ�����Сͣ��ʱ��

	double t_high = 400.0; //�´ε����������ʱ��
	double tf_low = 1.0 ;  //�����������С���ʱ��
	double p_factor = 1e4;

	//���ʱ���Ƿ��ͻ
	for (int i = 1; i < X.size(); i+=2)
	{
		tnow = X[i - 1];    //��һ����Ƭ�ĵ���ʱ��

		double penalt_temp = (X[i + 1] - tnow) - t_high;              //�����´ε����������ʱ�� tf_i+1 <tf_i + 30 
		if (penalt_temp > 0.0) penalty += p_factor * penalt_temp;
		
		double penalt_temp2 = t_low + tnow - X[i] ;                   //�����ٴγ�����Сͣ��ʱ�� tf_i + 5 <ts_i+1
		if (penalt_temp2 > 0.0) penalty += p_factor * penalt_temp2;

		double penalt_temp3 = X[i] + tf_low - X[i+1];                 //��������������С���ʱ�� ts_i+1 + 1 <tf_i+1
		if (penalt_temp3 > 0.0) penalty += p_factor * penalt_temp3;
	}

	double total_dv = Dv_All(X, node_order, derbrisnum);

	return total_dv + penalty;
}

/****************************************************************************
* ������   : map_input(vector<double>& X)
* ��  ��   : �����롾0,1��ӳ������ȷ��Χ
****************************************************************************/
void map_input(vector<double>& X)
{
	X[0] = X[0] * 500;
	for (int i= 1; i<X.size();i++)
	{
		X[i] = X[i] * 300;
		X[i] += X[i-1] ;
	}
}

/****************************************************************************
* ������   : nlopt_tstf_pso(const std::vector<double>& x, std::vector<double>& grad, void* f_data)
* ��  ��   : PSO�����õ��Ż��������ײ���������ʱ��2952������
****************************************************************************/
double nlopt_tstf_pso(const std::vector<double>& x, std::vector<double>& grad, void* f_data)
{

	vector<double> X(x);
	map_input(X);
	
	int* node_order = (int*)f_data;
	int derbrisnum = (X.size() + 1) / 2;

	double penalty = 0.0;
	double tnow = X[0]; //��ʼ����ʱ��
	double t_low = 5.0; //�ٴγ�����Сͣ��ʱ��

	double t_high = 800.0; //�´ε����������ʱ��
	double tf_low = 1.0;  //�����������С���ʱ��
	double p_factor = 1e4;

	//���ʱ���Ƿ��ͻ
	for (int i = 1; i < X.size(); i += 2)
	{
		tnow = X[i - 1];    //��һ����Ƭ�ĵ���ʱ��

		double penalt_temp = (X[i + 1] - tnow) - t_high;              //�����´ε����������ʱ�� tf_i+1 <tf_i + 30 
		if (penalt_temp > 0.0) penalty += p_factor * penalt_temp;

		double penalt_temp2 = t_low + tnow - X[i];                 //�����ٴγ�����Сͣ��ʱ�� tf_i + 5 <ts_i+1
		if (penalt_temp2 > 0.0) penalty += p_factor * penalt_temp2;

		double penalt_temp3 = X[i] + tf_low - X[i + 1];                 //��������������С���ʱ�� ts_i+1 + 1 <tf_i+1
		if (penalt_temp3 > 0.0) penalty += p_factor * penalt_temp3;
	}

	if( X.back()> 2952.0) penalty += p_factor *( X.back() - 2952.0);
	

	double total_dv = Dv_All(X, node_order, derbrisnum);

	return total_dv + penalty;
}

//mԼ��������nΪ�Ż���������
//Լ������С��0
double tstf_constrant(unsigned m, double* result, unsigned n, const double* X,double* gradient, void* func_data)
{

	return  0.0;
}