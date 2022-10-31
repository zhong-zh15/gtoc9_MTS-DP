/****************************************************************************
* Copyright (C), 2020-2031 Tsinghua University, School of Aerospace Engineering, LAD
* Author: Zhong Zhang
				zhong-zh19@mails.tsinghua.edu.cn
*               539977562@qq.com
* File: OrbitMath.h
* Description: basic math functions. This is developed by Fanghua Jiang
*
* Log:
*Version      Date        Author           Description
****************************************************************************/

#ifndef _ORBITMATH_H_
#define _ORBITMATH_H_
#include <math.h>
#include<assert.h>
#include"Constant.h"
#include <vector>

/****************************************************************************
* Function     : V_Copy
* Description  : assign the value of vector A to vector B
*                input:
*					A: vector
*					N: vector dimension
*                ouput:
*					B: vector 
****************************************************************************/
template<class T> inline void V_Copy(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++) B[I_]=A[I_];
}

/****************************************************************************
* Function     : V_Copy
* Description  : assign the three values x, y, z, to B in turn
*                input:
*					x, y, z: three values
*                ouput:
*					B: vector
****************************************************************************/
template<class T> inline void V_Copy(T* B, T x, T y, T z)
{
	B[0]=x;
	B[1]=y;
	B[2]=z;
}

/****************************************************************************
* Function     : V_Copy
* Description  : assign the six values x, y, z, vx, vy, vz to B in turn
*                input:
*					x, y, z, vx, vy, vz: six values
*                ouput:
*					B: vector
****************************************************************************/
template<class T> inline void V_Copy(T* B, T x, T y, T z, T vx, T vy, T vz)
{
	B[0]=x;
	B[1]=y;
	B[2]=z;
	B[3]=vx;
	B[4]=vy;
	B[5]=vz;
}

/****************************************************************************
* Function     : V_BoolEqua
* Description  : returns true if every element of vector B and vector A are equal
*                input:
*     				B: vector
*					A: vector 
*				    N: vector dimension
*                ouput:
*     				return value: 1 or 0
****************************************************************************/
template<class T> inline bool V_BoolEqua(const T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++)	if(B[I_]!=A[I_])return false;
	return true;
}

/****************************************************************************
* Function     : V_Opposite
* Description  : assign vector -A to vector B, B[i]=-A[i] 
*                input:
*     				A: vector
*					N: vector dimension
*                ouput:
*     				B: vector
****************************************************************************/
template<class T> inline void V_Opposite(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++) B[I_]=-A[I_];
}

/****************************************************************************
* Function     : V_Add
* Description  : assign vector A+B to C, C[i]=A[i]+B[i]
*                input:
*     				A: vector
*					B: vector
*					N: vector dimension
*                ouput:
*     				C: vector
****************************************************************************/
template<class T> inline void V_Add(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]+B[I_];	
}

/****************************************************************************
* Function     : V_Add
* Description  : each element of vector B plus number A equals vector C, C[i]=B[i]+A
*                input:
*					B: vector
*     				A: number
*					N: dimension of vector B
*                ouput:
*     				C: vector
****************************************************************************/
template<class T> inline void V_Add(T* C, const T* B, T A, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A+B[I_];	
}

/****************************************************************************
* Function     : V_Minus
* Description  : assign vector A-B to vector C, C[i]=A[i]-B[i]
*                input:
*     				A: vector
*					B: vector
*					N: vector dimension
*                ouput:
*     				C: vector
****************************************************************************/
template<class T> inline void V_Minus(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]-B[I_];	
}

/****************************************************************************
* Function     : V_Minus
* Description  : each element of vector A minus number B equals vector C, C[i]=A[i]-B
*                input:
*					A: vector
*     				B: number
*					N: dimension of vector A
*                ouput:
*     				C: vector
****************************************************************************/
template<class T> inline void V_Minus(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]-B;
}

/****************************************************************************
* Function     : V_Multi
* Description  : for vector C, C[i]=A[i]*B[i]
*                input:
*					A: vector
*     				B: vector
*					N: dimension of vector A, B 
*                ouput:
*     				C: vector
****************************************************************************/
template<class T> inline void V_Multi(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]*B[I_];	
}

/****************************************************************************
* Function     : V_Multi
* Description  : for vector C, C[i]=A[i]*B
*                input:
*					A: vector
*     				B: number
*					N: dimension of vector A
*                ouput:
*     				C: vector
****************************************************************************/
template<class T> inline void V_Multi(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=B*A[I_];
}

/****************************************************************************
* Function     : V_Divid
* Description  : for vector C, C[i]=A[i]/B[i]
*                input:
*					A: vector
*     				B: vector
*					N: dimension of vector A, B
*                ouput:
*     				C: vector
****************************************************************************/
template<class T> inline void V_Divid(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]/B[I_];	
}

/****************************************************************************
* Function     : V_Divid
* Description  : for vector C, C[i]=A[i]/B
*                input:
*					A: vector
*     				B: number
*					N: dimension of vector A
*                ouput:
*     				C: vector
****************************************************************************/
template<class T> inline void V_Divid(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]/B;
}

/****************************************************************************
* Function     : V_Dot
* Description  : find the inner product of vector A and vector B
*                input:
*					A: vector
*     				B: vector
*					N: dimension of vector A, B 
*                ouput:
*     				return value: value of inner product
****************************************************************************/
template<class T> inline T V_Dot(const T* A, const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;I_++) result+=A[I_]*B[I_];
	return result;
}

/****************************************************************************
* Function     : V_Cross
* Description  : find the outer product of vector A and vector B, cannot use V_Cross(B,B,A) or V_Cross(B,A,B)
*                input:
*					A: vector
*     				B: vector
*                ouput:
*     				C: vector, C=AXB
****************************************************************************/
template<class T> inline void V_Cross(T* C, const T* A, const T* B)
{
	C[0]=A[1]*B[2]-A[2]*B[1];
	C[1]=A[2]*B[0]-A[0]*B[2];
	C[2]=A[0]*B[1]-A[1]*B[0];
}

/****************************************************************************
* Function     : V_Absol
* Description  : find the absolute value of each element of vector B, 
*                input:
*     				B: vector
*					N: dimension of vector B
*                ouput:
*     				A: vector, A[i]=|B[i]|
****************************************************************************/
template<class T> inline void V_Absol(T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++)
	{
		if(B[I_]>=0)
			A[I_]=B[I_];
		else
			A[I_]=-B[I_];
	}
}

/****************************************************************************
* Function     : V_Norm1
* Description  : find the 1-norm of vector B,
*                input:
*     				B: vector
*					N: dimension of vector B
*                ouput:
*     				return value: value of 1-norm of vector B
****************************************************************************/
template<class T> inline T V_Norm1(const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;I_++)
	{
		if(B[I_]>=0)
			result+=B[I_];
		else
			result-=B[I_];
	}
	return result;
}

/****************************************************************************
* Function     : V_Norm2
* Description  : find the 2-norm of vector B,
*                input:
*     				B: vector
*					N: dimension of vector B
*                ouput:
*     				return value: value of 2-norm of vector B
****************************************************************************/
template<class T> inline T V_Norm2(const T* B, int N)
{
	T result=V_Dot(B,B,N);
	return sqrt(result);
}

/****************************************************************************
* Function     : V_NormInf
* Description  : find the infinite-norm of vector B,
*                input:
*     				B: vector
*					N: dimension of vector B
*                ouput:
*     				return value: value of infinite-norm of vector B
****************************************************************************/
template<class T> inline T V_NormInf(const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;I_++) 
	{
		if(B[I_]>=0) {if(B[I_]>result) result=B[I_];}
		else {if(-B[I_]>result) result=-B[I_];}
	}
	return result;
}

/****************************************************************************
* Function     : V_Max
* Description  : find the largest element in the vector B
*                input:
*     				B: vector
*					N: dimension of vector B
*                ouput:
*     				return value: value of the largest element
****************************************************************************/
template<class T> inline T V_Max(const T* B, int N)
{
	T result=B[0];
	for(int I_=0;I_<N;I_++) if(B[I_]>result) result=B[I_];
	return result;
}

/****************************************************************************
* Function     : V_Max
* Description  : find the maximum element in the vector B
*                input:
*     				B: vector
*					N: dimension of vector B
*                ouput:
*     				return value: value of the largest element
*				    index: numerical order of the largest element
****************************************************************************/
template<class T> inline T V_Max(int & index, const T* B, int N)
{
	T maximal=B[0];
	index=0;
	for(int I_=0;I_<N;I_++) if(B[I_]>maximal) {index=I_;maximal=B[I_];}
	return maximal;
}

/****************************************************************************
* Function     : V_Min
* Description  : find the minimum element in the vector B
*                input:
*     				B: vector
*					N: dimension of vector B
*                ouput:
*     				return value: value of the minimum element
****************************************************************************/
template<class T> inline T V_Min(const T* B, int N)
{
	T result=B[0];
	for(int I_=0;I_<N;I_++) if(B[I_]<result) result=B[I_];
	return result;
}

/****************************************************************************
* Function     : V_Min
* Description  : find the minimum element in the vector B
*                input:
*     				B: vector
*					N: dimension of vector B
*                ouput:
*     				return value: value of the minimum element
*				    index: numerical order of the minimum element
****************************************************************************/
template<class T> inline T V_Min(int& index, const T* B, int N)
{
	T minimal=B[0];
	index=0;
	for(int I_=0;I_<N;I_++) if(B[I_]<minimal) {index=I_;minimal=B[I_];}
	return minimal;
}

/****************************************************************************
* Function     : Sign
* Description  : sign function
*                input:
*     				InputValue: number
*                ouput:
*     				return value: 1 or -1 or 0
****************************************************************************/
template<class T> inline int Sign(const T & InputValue)
{
	if(InputValue>0)
		return 1;
	else if(InputValue<0)
		return -1;
	else
		return 0;
}

/****************************************************************************
* Function     : Max
* Description  : find the maximum of two numbers
*                input:
*     				x: number
*					y: number
*                ouput:
*     				return value: the maximum of x and y
****************************************************************************/
template <class T>
inline T Max (T x, T y) 
{ 
	return (x>y)?x:y;
}

/****************************************************************************
* Function     : Min
* Description  : find the minimum of two numbers
*                input:
*     				x: number
*					y: number
*                ouput:
*     				return value: the minimum of x and y
****************************************************************************/
template <class T>
inline T Min (T x, T y) 
{
	return (x<y)?x:y;
}

/****************************************************************************
* Function     : NiceAngle
* Description  : calculate the quadrant angle corresponding to the coordinates (x, y), [0,2*PI)
*                input:
*     				x: number
*					y: number
*                ouput:
*     				return value: the quadrant angle, [0,2*PI)
****************************************************************************/
inline double NiceAngle(double x, double y)
{
	double temp=atan2(y, x);
	if(temp<0.0)
		temp+=2*DPI;
	return temp;
}

/****************************************************************************
* Function     : NiceAngle
* Description  : convert any angle size to [0,2*PI)
*                input:
*     				alpha: angle (unit: rad)
*                ouput:
*     				return value: [0,2*PI)
****************************************************************************/
inline double NiceAngle(double alpha)
{
	double temp=fmod(alpha, 2.0*DPI);
	if(temp<0.0)
		temp+=2.0*DPI;
	return temp;
}

/****************************************************************************
* Function     : atanh
* Description  : atanh function
*                input:
*     				x: number
*                ouput:
*     				return value: atanh(x)
****************************************************************************/
inline double atanh(double x)
{
	assert(fabs(x)<1.0);
	return 0.5*log((1.0+x)/(1.0-x));
}

/****************************************************************************
* Function     : asinh
* Description  : asinh function
*                input:
*     				x: number
*                ouput:
*     				return value: asinh(x)
****************************************************************************/
inline double asinh(double x)
{
	
	return log(x+sqrt(x*x+1.0));
}

/****************************************************************************
* Function     : acosh
* Description  : acosh function
*                input:
*     				x: number
*                ouput:
*     				return value: acosh(x)
****************************************************************************/
inline double acosh(double x)
{	
	assert(x>=1.0);
	return log(x+sqrt(x*x-1.0));
}

#endif