#include "OrbitMath.h"


double LagInterp(double *x, double* y, int order, double t)
{
	if(order<=0)
		return y[0];
	if(order==1)
		return ((t-x[1])*y[0]-(t-x[0])*y[1])/(x[0]-x[1]);
	double temp=0;
	for(int k=0;k<=order;k++)
	{
		double temps=1.0;
		for(int i=0;i<=order;i++)
			if(i!=k) temps*=(t-x[i])/(x[k]-x[i]);
		temp+=temps*y[k];
	}
	return temp;
}