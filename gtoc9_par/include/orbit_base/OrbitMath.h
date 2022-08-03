#ifndef _ORBITMATH_H_
#define _ORBITMATH_H_
#include <math.h>
#include<assert.h>
#include"Constant.h"
#include <vector>

//将向量A的值赋给B
template<class T> inline void V_Copy(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++) B[I_]=A[I_];
}

//将三个值依次赋给B
template<class T> inline void V_Copy(T* B, T x, T y, T z)
{
	B[0]=x;
	B[1]=y;
	B[2]=z;
}

//将六个值依次赋给B
template<class T> inline void V_Copy(T* B, T x, T y, T z, T vx, T vy, T vz)
{
	B[0]=x;
	B[1]=y;
	B[2]=z;
	B[3]=vx;
	B[4]=vy;
	B[5]=vz;
}

//向量B与A每个元素都相等时返回真
template<class T> inline bool V_BoolEqua(const T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++)	if(B[I_]!=A[I_])return false;
	return true;
}

//B[i]=-A[i]
template<class T> inline void V_Opposite(T* B, const T* A, int N)
{
	for(int I_=0;I_<N;I_++) B[I_]=-A[I_];
}

//向量C[i]=A[i]+B[i]
template<class T> inline void V_Add(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]+B[I_];	
}

//向量C[i]=B[i]+A
template<class T> inline void V_Add(T* C, const T* B, T A, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A+B[I_];	
}

//向量C[i]=A[i]-B[i]
template<class T> inline void V_Minus(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]-B[I_];	
}

//向量C[i]=A[i]-B
template<class T> inline void V_Minus(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]-B;
}

//向量C[i]=A[i]*B[i]
template<class T> inline void V_Multi(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]*B[I_];	
}

//向量C[i]=A[i]*B
template<class T> inline void V_Multi(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=B*A[I_];
}

//向量C[i]=A[i]/B[i]
template<class T> inline void V_Divid(T* C, const T* A, const T* B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]/B[I_];	
}

//向量C[i]=A[i]/B
template<class T> inline void V_Divid(T* C, const T* A, T B, int N)
{
	for(int I_=0;I_<N;I_++) C[I_]=A[I_]/B;
}

//求内积
template<class T> inline T V_Dot(const T* A, const T* B, int N)
{
	T result=0;
	for(int I_=0;I_<N;I_++) result+=A[I_]*B[I_];
	return result;
}

//求外积C=AXB,不能用V_Cross(B,B,A)或V_Cross(B,A,B)
template<class T> inline void V_Cross(T* C, const T* A, const T* B)
{
	C[0]=A[1]*B[2]-A[2]*B[1];
	C[1]=A[2]*B[0]-A[0]*B[2];
	C[2]=A[0]*B[1]-A[1]*B[0];
}

//求A[i]=|B[i]|
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

//求1-范数
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

//求2-范数
template<class T> inline T V_Norm2(const T* B, int N)
{
	T result=V_Dot(B,B,N);
	return sqrt(result);
}

//求无穷-范数
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

//求最大元素
template<class T> inline T V_Max(const T* B, int N)
{
	T result=B[0];
	for(int I_=0;I_<N;I_++) if(B[I_]>result) result=B[I_];
	return result;
}

//求最大元素
template<class T> inline T V_Max(int & index, const T* B, int N)
{
	T maximal=B[0];
	index=0;
	for(int I_=0;I_<N;I_++) if(B[I_]>maximal) {index=I_;maximal=B[I_];}
	return maximal;
}

//求最小元素
template<class T> inline T V_Min(const T* B, int N)
{
	T result=B[0];
	for(int I_=0;I_<N;I_++) if(B[I_]<result) result=B[I_];
	return result;
}

//求最小元素
template<class T> inline T V_Min(int& index, const T* B, int N)
{
	T minimal=B[0];
	index=0;
	for(int I_=0;I_<N;I_++) if(B[I_]<minimal) {index=I_;minimal=B[I_];}
	return minimal;
}

//符号函数
template<class T> inline int Sign(const T & InputValue)
{
	if(InputValue>0)
		return 1;
	else if(InputValue<0)
		return -1;
	else
		return 0;
}

//求最大值
template <class T>
inline T Max (T x, T y) 
{ 
	return (x>y)?x:y;
}


//求最小值
template <class T>
inline T Min (T x, T y) 
{
	return (x<y)?x:y;
}

template <class T> inline T CopySign (T x, T y)
{
	return (y < 0) ? ((x < 0) ? x : -x) : ((x > 0) ? x : -x);
}

template <class T> void PaiXu(T* p, int N)
{

	T tmp = 0; 
	for (int i=0 ;i<N ;i++ ) 
	{ 
		for (int j=0 ;j<N-1-i ;j++ ) 
		{ 
			if (p[j] > p[j+1]) 
			{ 
				tmp = p[j]; 
				p[j] = p[j+1]; 
				p[j+1] = tmp;
			}
		} 
	}
}



//坐标(x,y)对应的相角,[0,2*PI).
inline double NiceAngle(double x, double y)
{
/*	if(x<=-EPSILON)
		return PI-atan(-y/x);
	else if(x>=EPSILON)
	{
		if(y<-EPSILON)
			return 2.0*PI+atan(y/x);
		else
			return atan(y/x);
	}	
	else
	{
		if(y>EPSILON)	return 0.5*PI;
		else if(y<-EPSILON) return 1.5*PI;
		else return 0.0;
	}*/
	double temp=atan2(y, x);
	if(temp<0.0)
		temp+=2*DPI;
	return temp;
}

//将角度值alpha转换到[0,2*PI)中.
inline double NiceAngle(double alpha)
{
/*	double x=cos(alpha);
	double y=sin(alpha);
	return NiceAngle(x, y);*/
	double temp=fmod(alpha, 2.0*DPI);
	if(temp<0.0)
		temp+=2.0*DPI;
	return temp;
}

inline double NiceAngle_pi2(double alpha)
{
	/*	double x=cos(alpha);
		double y=sin(alpha);
		return NiceAngle(x, y);*/
	double temp = fmod(alpha, 2.0 * DPI);
	if (temp < 0.0)
		temp += 2.0 * DPI;
	return temp;
}

//取最近接x的整数,返回double型.
inline double ANINT(double x)
{
	double left=fmod(x, 1.0);
	if(fabs(left)<0.5)
		return x-left;
	else if(left>=0.5)
		return x-left+1;
	else
		return x-left-1;
}

//取最近接x的整数,返回int型.
inline int NINT(double x)
{
	double left=ANINT(x);
	return (int)left;
}

inline double atanh(double x)
{
	assert(fabs(x)<1.0);
	return 0.5*log((1.0+x)/(1.0-x));
}

inline double asinh(double x)
{
	
	return log(x+sqrt(x*x+1.0));
}

inline double acosh(double x)
{	
	assert(x>=1.0);
	return log(x+sqrt(x*x-1.0));
}

//二维指针数组表示
//释放内存
template<class T> inline void M_Dele(T** B, int N, int M)
{
	for(int I_=0;I_<N;I_++) {delete[] B[I_];B[I_]=NULL;}
	delete[] B;
	B=NULL;
}
//将A的值赋给B
template<class T> inline void M_Copy(T** B, T** A, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) B[I_][J_]=A[I_][J_];
}
//将九个值依次赋给B
template<class T> inline void M_Copy(T** B, T a11, T a12, T a13, T a21, T a22, T a23, T a31, T a32, T a33)
{
	B[0][0]=a11;
	B[0][1]=a12;
	B[0][2]=a13;
	B[1][0]=a21;
	B[1][1]=a22;
	B[1][2]=a23;
	B[2][0]=a31;
	B[2][1]=a32;
	B[2][2]=a33;
}

//将数组A的值逐个赋给B
template<class T> inline void M_Copy(T** B, T* A, int N, int M)
{
	int i=0;
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) B[I_][J_]=A[i++];
}
//矩阵B与A每个元素都相等时返回真
template<class T> inline bool M_BoolEqua(T** B, T** A, int N, int M)
{
	for(int I_=0;I_<N;I_++)	for(int J_=0;J_<M;J_++) if(B[I_][J_]!=A[I_][J_])return false;
	return true;
}
//B[i][j]=-A[i][j]
template<class T> inline void M_Opposite(T** B, T** A, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) B[I_][J_]=-A[I_][J_];
}
//B[i][j]=A[j][i]
template<class T> inline void M_Tranpose(T** B, T** A, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) B[I_][J_]=A[J_][I_];
}

//C[i][j]=A[i][j]+B[i][j]
template<class T> inline void M_Add(T** C, T** A, T** B, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) C[I_][J_]=A[I_][J_]+B[I_][J_];	
}

//C[i][j]=B[i][j]+A
template<class T> inline void M_Add(T** C, T** B, T A, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) C[I_][J_]=A+B[I_][J_];	
}
//C[i][j]=A[i][j]-B[i][j]
template<class T> inline void M_Minus(T** C, T** A, T** B, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) C[I_][J_]=A[I_][J_]-B[I_][J_];	
}
//C[i][j]=A[i][k]*B[k][j]
template<class T> inline void M_Multi(T** C, T** A, T** B, int N, int M, int K)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++)
	{
		C[I_][J_]=0;
		for(int i=0;i<K;i++){ C[I_][J_]+=A[I_][i]*B[i][J_];	}
	}
}
//C[i]=A[i][j]*B[j]
template<class T> inline void M_Multi(T* C, T** A, T* B, int N, int M)
{
	for(int I_=0;I_<N;I_++)
	{
		C[I_]=0;
		for(int J_=0;J_<M;J_++)	C[I_]+=A[I_][J_]*B[J_];
	}
}
//C[i][j]=A[i][j]*B
template<class T> inline void M_Multi(T** C, T** A, T B, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) C[I_][J_]=B*A[I_][J_];
}


//C[i][j]=A[i][j]/B[i][j]
template<class T> inline void M_Divid(T** C, T** A, T** B, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) C[I_][J_]=A[I_][J_]/B[I_][J_];	
}

//求A[i][j]=|B[i][j]|
template<class T> inline void M_Absol(T** A, T** B, int N, int M)
{
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++)
	{
		if(B[I_][J_]>=0)
			A[I_][J_]=B[I_][J_];
		else
			A[I_][J_]=-B[I_][J_];
	}
}
//求最大元素
template<class T> inline T M_Max(T** B, int N, int M)
{
	T result=B[0][0];
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) if(B[I_][J_]>result) result=B[I_][J_];
	return result;
}


//求最小元素
template<class T> inline T M_Min(T** B, int N, int M)
{
	T result=B[0][0];
	for(int I_=0;I_<N;I_++) for(int J_=0;J_<M;J_++) if(B[I_][J_]<result) result=B[I_][J_];
	return result;
}

//求逆
template<class T> inline void swaprows(T** A, int row0, int row1)
{
    T* temp=NULL;
    temp=A[row0];
    A[row0]=A[row1];
    A[row1]=temp;
}

//	gjelim 
void M_Inverse(double** B, double** A, int dim, double** wa); 

////Lagrange插值函数
double LagInterp(double *x, double* y, int order, double t);

//二分查找
template <typename T> inline int binSearch(const std::vector<T> &S, T const& e, int lo, int hi)
{
	while (lo < hi)
	{
		int mi = (lo + hi) >> 1;
		e < S[mi] ? hi = mi : lo = mi + 1;
	}
	return --lo;
}
//二分查找
template <typename T> inline int binSearch(T* S, T const& e, int lo, int hi)
{
	while (lo < hi)
	{
		int mi = (lo + hi) >> 1;
		e < S[mi] ? hi = mi : lo = mi + 1;
	}
	return --lo;
}

#endif