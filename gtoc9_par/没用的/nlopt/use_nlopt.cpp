#include "use_nlopt.h"


/****************************************************************************
* 函数名   : nloptmain()
* 功  能   : 使用nlopt优化ts，tf，同时需要考虑约束（使用罚函数）
****************************************************************************/
void nloptmain(double (*ObjFun)(const std::vector<double>& x, std::vector<double>& grad, void* f_data), int* node_order,
	std::vector<double>& xbest, double& fbest, int ItMax)
{
	int num_variable = xbest.size();                                //只优化tf

	nlopt::opt opter(nlopt::LN_SBPLX, num_variable);				//局部优化  LN_SBPLX  LN_COBYLA LN_BOBYQA
		                                                                      //GN_DIRECT x  GD_STOGO x  GN_AGSx GN_CRS2_LM GN_ISRESx GN_ESCH
		                                                                      //效果来看 LN_SBPLX、 GN_CRS2_LM  比较好

	opter.set_min_objective(ObjFun, node_order);					//指标
	
	std::vector<double> lb(num_variable), rb(num_variable), dx(num_variable);

	//opter.add_inequality_mconstraint(ObjFun, node_order, dx);
	
	fbest = 1.0e20;
	for (int i = 0; i < num_variable; i++)
	{
		lb[i] = 0.0;												//下界
		rb[i] = 2952.0;												//上界
		dx[i] = 1.0/86400.0*200.0;									//初始步长
	}
	//opter.set_initial_step(dx);										//设置初始步长
	opter.set_lower_bounds(lb);										//设置下届
	opter.set_upper_bounds(rb);										//设置上届
	double tol = 1e-8;
	//opter.set_ftol_abs(tol);
	//opter.set_force_stop(tol);

	opter.set_maxeval(ItMax);										//优化该次数后停止

	nlopt::result res = opter.optimize(xbest, fbest);				//进行优化

}


/****************************************************************************
* 函数名   : nlopt_tf()
* 功  能   : 使用nlopt优化tf(n-1个优化变量，ts直接生成)，同时需要考虑约束（使用罚函数）上界和下届
****************************************************************************/
double nlopt_tstf(const std::vector<double>& X, std::vector<double>& grad, void* f_data)
{
	int* node_order = (int*)f_data;
	int derbrisnum = (X.size() + 1 )/2 ;

	double penalty = 0.0;
	double tnow =  X[0]; //初始出发时刻
	double t_low = 5.0; //再次出发最小停留时间

	double t_high = 400.0; //下次到达最大任务时间
	double tf_low = 1.0 ;  //出发到达的最小间距时间
	double p_factor = 1e4;

	//检测时间是否冲突
	for (int i = 1; i < X.size(); i+=2)
	{
		tnow = X[i - 1];    //上一个碎片的到达时刻

		double penalt_temp = (X[i + 1] - tnow) - t_high;              //满足下次到达最大任务时间 tf_i+1 <tf_i + 30 
		if (penalt_temp > 0.0) penalty += p_factor * penalt_temp;
		
		double penalt_temp2 = t_low + tnow - X[i] ;                   //满足再次出发最小停留时间 tf_i + 5 <ts_i+1
		if (penalt_temp2 > 0.0) penalty += p_factor * penalt_temp2;

		double penalt_temp3 = X[i] + tf_low - X[i+1];                 //满足出发到达的最小间距时间 ts_i+1 + 1 <tf_i+1
		if (penalt_temp3 > 0.0) penalty += p_factor * penalt_temp3;
	}

	double total_dv = Dv_All(X, node_order, derbrisnum);

	return total_dv + penalty;
}

/****************************************************************************
* 函数名   : map_input(vector<double>& X)
* 功  能   : 将输入【0,1】映射至正确范围
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
* 函数名   : nlopt_tstf_pso(const std::vector<double>& x, std::vector<double>& grad, void* f_data)
* 功  能   : PSO所调用的优化函数，底部加入了总时长2952的限制
****************************************************************************/
double nlopt_tstf_pso(const std::vector<double>& x, std::vector<double>& grad, void* f_data)
{

	vector<double> X(x);
	map_input(X);
	
	int* node_order = (int*)f_data;
	int derbrisnum = (X.size() + 1) / 2;

	double penalty = 0.0;
	double tnow = X[0]; //初始出发时刻
	double t_low = 5.0; //再次出发最小停留时间

	double t_high = 800.0; //下次到达最大任务时间
	double tf_low = 1.0;  //出发到达的最小间距时间
	double p_factor = 1e4;

	//检测时间是否冲突
	for (int i = 1; i < X.size(); i += 2)
	{
		tnow = X[i - 1];    //上一个碎片的到达时刻

		double penalt_temp = (X[i + 1] - tnow) - t_high;              //满足下次到达最大任务时间 tf_i+1 <tf_i + 30 
		if (penalt_temp > 0.0) penalty += p_factor * penalt_temp;

		double penalt_temp2 = t_low + tnow - X[i];                 //满足再次出发最小停留时间 tf_i + 5 <ts_i+1
		if (penalt_temp2 > 0.0) penalty += p_factor * penalt_temp2;

		double penalt_temp3 = X[i] + tf_low - X[i + 1];                 //满足出发到达的最小间距时间 ts_i+1 + 1 <tf_i+1
		if (penalt_temp3 > 0.0) penalty += p_factor * penalt_temp3;
	}

	if( X.back()> 2952.0) penalty += p_factor *( X.back() - 2952.0);
	

	double total_dv = Dv_All(X, node_order, derbrisnum);

	return total_dv + penalty;
}

//m约束个数，n为优化变量个数
//约束满足小于0
double tstf_constrant(unsigned m, double* result, unsigned n, const double* X,double* gradient, void* func_data)
{

	return  0.0;
}