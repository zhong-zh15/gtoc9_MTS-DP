#include "OrbitFun.h"

#include <math.h>       
#include <stdio.h>


#include "Constant.h"
#include "OrbitMath.h"


//���������� coe[5]:
// 	coe[0]�볤��a���ף�:������Բ(��Բ)��˫���߹��,��ֵ������;���������߹��,
//                      ��ֵȡΪ���Ǿ�.
// 	coe[1]ƫ����e:����Բ���e=0;��Բ���0<e<1,�����߹��e=1,˫���߹��e>1
// 	coe[2]������i�����ȣ�:��Χ0<=i<=180��.
//	coe[3]������ྭOmega�����ȣ�:��������Ϊ0ʱû������,�ɽ���ֵȡΪ0
//	coe[4]���������omega�����ȣ�:��ƫ����Ϊ0ʱû������,�ɽ���ֵ��Ϊ0
//	coe[5]������f�����ȣ�:��ƫ����Ϊ0ʱû������,�ɽ���ֵȡΪomega+f,��γ�ȷ���

/*********************************�������������ֽǶȹ�ϵ***************************************************************/
// E2f ����ƫ����Ǻ�ƫ������������
//�������E, e:
//   E:ƫ�����,��λ������,����˫�����,ָH,��r=a(ecoshH-1)
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//�������f:
// 	�����ǣ����ȣ�
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
double E2f(int& flag, double E, double e)
{
	if (e < 0.0) { flag = 0; return E; }

	double f = 0.0;
	if (e >= 0.0 && e < 1.0)//Բ����Բ���
	{
		double E0 = fmod(E, D2PI);
		if (E0 > DPI)
			E0 -= D2PI;
		if (E0 < -DPI)
			E0 += D2PI;
		f = 2.0 * atan(sqrt((1.0 + e) / (1.0 - e)) * tan(0.5 * E0));
		f = f + E - E0;
	}
	else if (e > 1.0)//˫���߹��
		f = 2.0 * atan(sqrt((e + 1.0) / (e - 1.0)) * tanh(0.5 * E));
	else// (abs(e-1.0)<epsilon)�����߹��
	{
		f = E;
		//		cout<<"�����߹��û�ж���ƫ�����.�ڴ˽������ǵ�ֵ��ΪE."<<endl;
	}
	flag = 1;
	return f;
}
// E2M ����ƫ����Ǻ�ƫ������ƽ�����
//�������E, e:
//   E:ƫ�����,��λ������,����˫�����,ָH
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
//�������M:
// 	ƽ�����(����),����˫�����,ָN
double E2M(int& flag, double E, double e)
{
	if (e < 0.0) { flag = 0; return E; }
	double M = 0.0;
	if (e >= 0.0 && e < 1.0)//Բ����Բ���
	{
		double E0 = fmod(E, D2PI);
		M = E0 - e * sin(E0);
		M = M + E - E0;
	}
	else if (e > 1.0)//˫���߹��
		M = e * sinh(E) - E;
	else//(abs(e-1.0)<epsilon)�����߹��
	{
		M = E;
		//		cout<<"�����߹��û�ж���ƫ�����.�ڴ˽�ƽ�����ֵ��ΪE."<<endl;
	}
	flag = 1;
	return M;
}
// f2E ���������Ǻ�ƫ������ƫ�����
//�������f, e:
//   f:������,��λ������
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
//�������E:
// 	ƫ����ǣ����ȣ�,����˫�����,ָH,��r=a(ecoshH-1)
double f2E(int& flag, double f, double e)
{
	if (e < 0.0) { flag = 0; return f; }

	double E = 0.0;
	if (e >= 0.0 && e < 1.0)//Բ����Բ���
	{
		double f0 = fmod(f, D2PI);
		if (f0 > DPI)
			f0 -= D2PI;
		if (f0 < -DPI)
			f0 += D2PI;
		E = 2.0 * atan(sqrt((1.0 - e) / (1.0 + e)) * tan(0.5 * f0));
		E += f - f0;
	}
	else if (e > 1.0)//˫���߹��
	{
		if (f > DPI - acos(1.0 / e) || f < -DPI + acos(1.0 / e))
		{
			//			cout<<"�����ܴﵽ��˫�����."<<endl;
			flag = 0;
			return f;
		}
		else
			E = 2.0 * atanh(sqrt((e - 1.0) / (1.0 + e)) * tan(0.5 * f));
	}
	else//if(abs(e-1.0)<epsilon)�����߹��
	{
		E = f;
		//		cout<<"�����߹��û�ж���ƫ�����.�ڴ˽���ֵ��Ϊf."<<endl;
	}
	flag = 1;
	return E;
}
// M2E ����ƽ����Ǻ�ƫ������ƫ�����
//�������M, e, MaxIter, epsilon:
//   M:ƽ�����,��λ������,����˫�����,ָN
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//   MaxIter:����������,Ĭ��Ϊ60
//   epsilon:�����������,Ĭ��Ϊ1e-14;
//�������E:
// 	ƫ����ǣ����ȣ�����˫�������ָH
double M2E(int& flag, double M, double e, int MaxIter, double epsilon)
{
	if (epsilon <= 0.0 || MaxIter < 1 || e < 0.0) { flag = 0; return M; }

	//������������Solar System Dynamics��Chapter2,Carl D.Murray and Stanley F.Dermott��
	double E = 0.0, Minus = 0.0, DeMinus = 0.0, DeDeMinus = 0.0, DeDeDeMinus = 0.0, Delta1 = 0.0, Delta2 = 0.0, Delta3 = 0.0;
	int N = 0;
	if (e >= 0.0 && e < 1.0)//Բ����Բ���
	{
		double RM = fmod(M, D2PI);
		if (RM < 0.0)
			RM += D2PI;
		double sinRM = sin(RM);
		E = RM + 0.85 * e * Sign(sinRM);
		N = 0;
		Delta3 = 1.0;
		while (fabs(Delta3) >= epsilon && N < MaxIter)
		{
			Minus = E - e * sin(E) - RM;
			DeMinus = 1.0 - e * cos(E);
			DeDeMinus = e * sin(E);
			DeDeDeMinus = e * cos(E);
			Delta1 = -Minus / DeMinus;
			Delta2 = -Minus / (DeMinus + 0.5 * Delta1 * DeDeMinus);
			Delta3 = -Minus / (DeMinus + 0.5 * Delta2 * DeDeMinus + 1.0 / 6.0 * Delta2 * Delta2 * DeDeDeMinus);
			E = E + Delta3;
			N = N + 1;
		}
		E = E + M - RM;
	}
	else if (e > 1.0)//˫���߹��
	{
		E = asinh(M / e);
		Delta3 = 1.0;
		N = 0;
		while (fabs(Delta3) >= epsilon && N < MaxIter)
		{
			Minus = e * sinh(E) - E - M;
			DeMinus = e * cosh(E) - 1.0;
			DeDeMinus = e * sinh(E);
			DeDeDeMinus = e * cosh(E);
			Delta1 = -Minus / DeMinus;
			Delta2 = -Minus / (DeMinus + 0.5 * Delta1 * DeDeMinus);
			Delta3 = -Minus / (DeMinus + 0.5 * Delta2 * DeDeMinus + 1.0 / 6.0 * Delta2 * Delta2 * DeDeDeMinus);
			E = E + Delta3;
			N = N + 1;
		}
	}
	else //(abs(e-1.0)<epsilon)�����߹��
	{
		E = M;
		//		cout<<"�����߹��û�ж���ƫ�����.�ڴ˽���ֵ��ΪM."<<endl;
	}
	if (((e >= 0.0 && e < 1.0) || (e > 1.0)) && fabs(Delta3) >= 5.0 * epsilon && N >= MaxIter)
	{
		//		cout<<"����������,�뽵�;���epsilon�����ӵ�����������."<<endl;
		flag = 0;
		return M;
	}
	flag = 1;
	return E;
}
// f0dt2ft ���ݳ�ʼ�����Ǻ��ݻ�ʱ��������������
//�������f0,t,a,e,mu,MaxIter,epsilon:
//   f0:��ʼ������,��λ������
//   t:����ʱ��,��λ����
//   a:����볤�ᣬ��λ���ס����������߹����ָ���Ǿ࣬��p/2
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14
//   MaxIter:����������
//   epsilon:�����������,Ĭ��Ϊ1e-14;
//�������ft:
// 	����������(����)
double f0dt2ft(int& flag, double f0, double dt, double a, double e, double mu, int MaxIter, double epsilon)
{
	if (mu <= 0.0 || MaxIter < 1 || a <= 0.0 || e < 0.0) { flag = 0; return f0; }

	double ft = 0.0;
	if ((e >= 0.0 && e < 1.0) || (e > 1.0))//Բ,��Բ,˫�����
	{
		double E = f2E(flag, f0, e);
		if (flag == 0) return f0;
		double M = E2M(flag, E, e);
		if (flag == 0) return f0;
		M += sqrt(mu / (a * a * a)) * dt;
		E = M2E(flag, M, e, MaxIter, epsilon);
		if (flag == 0) return f0;
		ft = E2f(flag, E, e);
		if (flag == 0) return f0;
	}
	else //(abs(e-1.0)<epsilon)�����߹��
	{
		if ((f0 < -DPI) || (f0 > DPI))
		{
			//			cout<<"���������߹������ʼ������Ӧ��-180��180��֮��."<<endl;
			flag = 0; return f0;
		}
		else if (f0 > DPI || f0 < -DPI)
			ft = f0;
		else
		{
			double B = 0.75 * sqrt(2.0 * mu / (a * a * a)) * dt + 0.5 * tan(0.5 * f0) * ((tan(0.5 * f0)) * (tan(0.5 * f0)) + 3.0);
			double B1B = B + sqrt(1.0 + B * B);
			double tanv = 0.0;
			if (fabs(dt) < D2PI * sqrt((a * a * a) / mu) / 1000.0)//�ƽ�ʱ��ΪС�������
			{
				double A = pow(B1B, 2.0 / 3.0);
				tanv = 2.0 * A * B / (1.0 + (1.0 + A) * A);
			}
			else//����С�������
			{
				double temp = pow(B1B, 1.0 / 3.0);
				tanv = temp - 1.0 / temp;
			}
			ft = 2.0 * atan(tanv);
		}
	}
	flag = 1;
	return ft;
}
// f0ft2dt ���ݳ�ʼ�����Ǻ��������������ݻ�ʱ��
//�������f0,ft,a,e,mu,AbsTol:
//   f0:��ʼ������,��λ������
//   ft:�ն�������,��λ������
//   a:����볤�ᣬ��λ���ס����������߹����ָ���Ǿ࣬��p/2
//   e:ƫ����,0Բ,(0,1)��Բ,1������,(1,...)˫����
//   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
//      Ĭ��Ϊ��������ϵ��3.98600441800e+14
//   epsilon:��С����,Ĭ��Ϊ1.0e-12
//�������t:
//   t:����ƽ�ʱ��,��λ����
double f0ft2dt(int& flag, double f0, double ft, double a, double e, double mu)
{
	flag = 0;
	if (mu <= 0.0 || a <= 0.0 || e < 0.0 ) //|| ft < f0)
		return 0.0;

	double dt = 0.0;
	if (e >= 1.0)
	{
		double maxangle = DPI - acos(1.0 / e);
		if ((Min(f0, ft) < -maxangle) || (Max(f0, ft) > maxangle))
		{
			//			cout<<"�����ܴﵽ��˫���������߹��"<<endl;			
			return 0.0;
		}
		else if (f0<-maxangle || f0>maxangle || ft<-maxangle || ft>maxangle)
		{
			//			cout<<"����ʱ����Ѽ���׼ȷ,����Ϊ����,�ڴ���Ϊ1.0e100.";
			dt = 1.0e308;
			return dt;
		}
	}

	double omega = sqrt(mu / (a * a * a));
	double delta = 0.0;
	if ((e >= 0.0 && e < 1.0) || (e > 1.0))
	{
		double E = f2E(flag, f0, e);
		double M0 = E2M(flag, E, e);
		E = f2E(flag, ft, e);
		double Mt = E2M(flag, E, e);
		if (flag == 0) return 0.0;
		delta = Mt - M0;
	}
	else// if(fabs(e-1.0)<epsilon)
	{
		double B1 = tan(0.5 * f0) * ((tan(0.5 * f0)) * (tan(0.5 * f0)) + 3.0);
		double B2 = tan(0.5 * ft) * ((tan(0.5 * ft)) * (tan(0.5 * ft)) + 3.0);
		delta = sqrt(2.0) / 3.0 * (B2 - B1);
	}
	dt = delta / omega;
	flag = 1;
	return dt;
}

/**************************************************************************************************************************************/
/*******************************************************���ֵ��������Ƕȹ�ϵ*********************************************************/
/**************************************************************************************************************************************/
double L0dt2Lt(int& flag, double L0, double dt, const double* ee, double mu)
{
	double e = sqrt(ee[1] * ee[1] + ee[2] * ee[2]);
	double a = ee[0] / fabs(1.0 - e * e);
	double Oo = atan2(ee[2], ee[1]);
	double Lt = Oo + f0dt2ft(flag, L0 - Oo, dt, a, e, mu);
	return Lt;
}
double L0Lt2dt(int& flag, double L0, double Lt, const double* ee, double mu)
{
	double e = sqrt(ee[1] * ee[1] + ee[2] * ee[2]);
	double a = ee[0] / fabs(1.0 - e * e);
	double Oo = atan2(ee[2], ee[1]);
	double dt = f0ft2dt(flag, L0 - Oo, Lt - Oo, a, e, mu);
	return dt;
}

/*********************************������������ֱ�����ꡢ�Ľ����ֵ���������ת��**************************************/
void coe2rv(int& flag, double* rv, const double* coe, double mu = 3.98600441800e+14)
{
	flag = 0;
	//if (mu <= 0.0 || coe[0] <= 0.0 || coe[1] < 0.0 || coe[2]<0.0 || coe[2]>DPI)
	//	return;
	//if ((coe[1] * cos(coe[5])) < -1.0)
	//{
	//	//		cout<<"�����ܴﵽ��˫�����."<<endl;		
	//	return;
	//}

	double p = coe[0] * fabs(1.0 - coe[1] * coe[1]);//��ͨ��
	//if (coe[1] == 1.0)//����������߹��,����Դ�.
	//	p = 2.0 * coe[0];

	double sini, cosi, sinO, cosO, sino, coso;
	sini = sin(coe[2]);
	cosi = cos(coe[2]);
	sinO = sin(coe[3]);
	cosO = cos(coe[3]);
	sino = sin(coe[4]);
	coso = cos(coe[4]);

	//���ƽ�淨��λʸ��,���Ƕ�����λʸ��
	double HVector[3] = { sini * sinO, -sini * cosO, cosi };

	//ƫ���ʵ�λʸ��,���Laplaceʸ��
	double PVector[3] = { cosO * coso - sinO * sino * cosi, sinO * coso + cosO * sino * cosi, sino * sini };

	//��ͨ������λʸ��,PVector,QVector,HVector������������ϵ
	// QVector=[-cosO*sino-sinO*coso*cosi;-sinO*sino+cosO*coso*cosi;coso*sini];
	double QVector[3];
	V_Cross(QVector, HVector, PVector);

	double r = 0.0;
	//if ((coe[1] * cos(coe[5])) + 1.0 <= 0.0)
	//{
	//	//		cout<<"�����˫�������������Զ��."<<endl;
	//	r = 1.0e308;
	//}
	//else
		r = p / (1.0 + coe[1] * cos(coe[5]));

	for (int i = 0; i < 3; i++)
	{
		rv[i] = r * (cos(coe[5]) * PVector[i] + sin(coe[5]) * QVector[i]);
		rv[3 + i] = sqrt(mu / p) * (-sin(coe[5]) * PVector[i] + (cos(coe[5]) + coe[1]) * QVector[i]);
	}
	flag = 1;
	return;
}
void rv2coe(int& flag, double* coe, const double* RV, double mu = 3.98600441500e+14)
{
	int i;
	flag = 0;
	//if (mu <= 0.0)
	//	return;

	double R[3] = { RV[0], RV[1], RV[2] };
	double V[3] = { RV[3], RV[4], RV[5] };
	double radius = V_Norm2(R, 3);//����
	double velocity = V_Norm2(V, 3);//�ٶ�
	//if (radius <= 0.0 || velocity <= 0.0)
	//	return;
	double unitR[3];
	for (i = 0; i < 3; i++) unitR[i] = R[i] / radius;//����λʸ��    
	double unitV[3];
	for (i = 0; i < 3; i++) unitV[i] = V[i] / velocity;//����λʸ��
	double hvector[3];
	V_Cross(hvector, unitR, unitV);
	double h = radius * velocity * V_Norm2(hvector, 3);//�Ƕ���ֵ
	//if (h <= 0.0)
	//	return;
	double unith[3];
	for (i = 0; i < 3; i++) unith[i] = hvector[i] / V_Norm2(hvector, 3);//����淨��λʸ��
	//ƫ����ʸ��
	double evector[3];
	V_Cross(evector, unitV, unith);
	for (i = 0; i < 3; i++) evector[i] = (velocity * h / mu) * evector[i] - unitR[i];
	coe[1] = V_Norm2(evector, 3);//ƫ����
	double p = h * h / mu;
	//if (coe[1] == 1.0)
	//	coe[0] = 0.5 * p;//�����߹���Ľ��Ǿ�
	//else
		coe[0] = p / (fabs(1.0 - coe[1] * coe[1]));//�볤��
	bool judge = (coe[1] > 0.0);
	double unite[3] = { 0.0 };
	if (judge)
		for (i = 0; i < 3; i++) unite[i] = evector[i] / coe[1];//ƫ���ʵ�λʸ��
	coe[2] = acos(unith[2]);//������

	double unitN[3] = { -unith[1], unith[0], 0.0 };//����ʸ��,δ��һ��

	double temp[3];

	if (V_Norm2(unitN, 3) == 0.0)
	{
		coe[3] = 0.0;//������ྭ
//		cout<<"�����ǽӽ�0��180��,������ྭ��������.�ڴ˽�����Ϊ��."<<endl;
		if (!judge)
		{
			coe[4] = 0.0;//���ǵ����
//			cout<<"ƫ���ʽӽ�0,���ǵ������������.�ڴ˽�����Ϊ��."<<endl;        
			coe[5] = atan2(unitR[1] * unith[2], unitR[0]);//������
		}
		else
		{
			V_Cross(temp, unite, unitR);
			coe[4] = atan2(unite[1] * unith[2], unite[0]); //���ǵ����       
			coe[5] = atan2(V_Dot(unith, temp, 3), V_Dot(unite, unitR, 3));
		}
	}
	else
	{
		V_Cross(temp, unitN, unitR);
		coe[3] = atan2(unith[0], -unith[1]);
		coe[5] = atan2(V_Dot(unith, temp, 3), V_Dot(unitN, unitR, 3));
		if (!judge)
		{
			coe[4] = 0.0;
			//			cout<<"ƫ���ʽӽ�0,���ǵ������������.�ڴ˽�����Ϊ��."<<endl;
		}
		else
		{
			V_Cross(temp, unitN, unite);
			coe[4] = atan2(V_Dot(unith, temp, 3), V_Dot(unite, unitN, 3));
			coe[5] = coe[5] - coe[4];
		}
	}
	//ת����[0,2pi)��
	coe[3] = fmod(coe[3], D2PI);
	if (coe[3] < 0.0)
		coe[3] += D2PI;
	coe[4] = fmod(coe[4], D2PI);
	if (coe[4] < 0.0)
		coe[4] += D2PI;
	coe[5] = fmod(coe[5], D2PI);
	if (coe[1] >= 1.0)
	{
		if (coe[5] > DPI - acos(1.0 / coe[1]))
			coe[5] -= D2PI;
		else if (coe[5] < -DPI + acos(1.0 / coe[1]))
			coe[5] += D2PI;
	}
	flag = 1;
	return;
}
void coe2ee(int& flag, double* ee, const double* coe, double mu = 3.98600441500e+14)
{
	flag = 0;
	if (mu <= 0.0 || coe[0] < 0.0 || coe[1] < 0.0 || coe[2]<0.0 || coe[2]>DPI)
		return;
	if ((coe[1] * cos(coe[5])) < -1.0)
		//		cout<<"�����ܴﵽ��˫�����."<<endl;
		return;

	ee[0] = coe[0] * fabs(1.0 - coe[1] * coe[1]);//��ͨ��p
	if (coe[1] == 1.0)//����������߹��,����Դ�.
		ee[0] = 2.0 * coe[0];
	ee[1] = coe[1] * cos(coe[4] + coe[3]);//f
	ee[2] = coe[1] * sin(coe[4] + coe[3]);//g
	double temp = tan(coe[2] / 2.0);
	ee[3] = temp * cos(coe[3]);//h
	ee[4] = temp * sin(coe[3]);//k
	ee[5] = coe[4] + coe[3] + coe[5];
	ee[5] = fmod(ee[5], D2PI);
	if (ee[5] < 0.0)
		ee[5] += D2PI;
	flag = 1;
	return;
}
void ee2coe(int& flag, double* coe, const double* ee, double mu = 3.98600441800e+14)
{
	flag = 0;
	if (mu <= 0.0 || ee[0] <= 0.0)
		return;

	double p = ee[0], f = ee[1], g = ee[2], h = ee[3], k = ee[4], L = ee[5];
	coe[1] = sqrt(f * f + g * g);
	if (coe[1] == 1.0)
		coe[0] = 0.5 * p;//�����߹���Ľ��Ǿ�
	else
		coe[0] = p / (fabs(1.0 - coe[1] * coe[1]));//�볤��
	double temp = sqrt(h * h + k * k);
	coe[2] = 2.0 * atan(temp);
	if (temp <= 0.0)
	{
		coe[3] = 0.0;//������ྭ
//		cout<<"�����ǽӽ�0��180��,������ྭ��������.�ڴ˽�����Ϊ��."<<endl;
		if (coe[1] <= 0.0)
		{
			coe[4] = 0.0;//���ǵ����
//			cout<<"ƫ���ʽӽ�0,���ǵ������������.�ڴ˽�����Ϊ��."<<endl;        
			coe[5] = L;//������
		}
		else
		{
			coe[4] = atan2(g, f); //���ǵ����       
			coe[5] = L - coe[4];
		}
	}
	else
	{
		coe[3] = atan2(k, h);
		coe[5] = L - coe[3];
		if (coe[1] <= 0.0)
		{
			coe[4] = 0.0;
			//			cout<<"ƫ���ʽӽ�0,���ǵ������������.�ڴ˽�����Ϊ��."<<endl;
		}
		else
		{
			coe[4] = atan2(g * h - f * k, f * h + g * k);
			coe[5] = coe[5] - coe[4];
		}
	}
	//ת����[0,2pi)��
	coe[3] = fmod(coe[3], D2PI);
	if (coe[3] < 0.0)
		coe[3] += D2PI;
	coe[4] = fmod(coe[4], D2PI);
	if (coe[4] < 0.0)
		coe[4] += D2PI;
	coe[5] = fmod(coe[5], D2PI);
	if (coe[5] < 0.0)
		coe[5] += D2PI;
	if (coe[1] >= 1.0)
	{
		if (coe[5] > DPI - acos(1.0 / coe[1]))
			coe[5] -= D2PI;
		else if (coe[5] < -DPI + acos(1.0 / coe[1]))
			coe[5] += D2PI;
	}
	flag = 1;
	return;
}
void ee2rv(int& flag, double* rv, const double* ee, double mu = 3.98600441800e+14)
{
	flag = 0;
	if (mu <= 0.0 || ee[0] <= 0.0)
		return;

	double p = ee[0], f = ee[1], g = ee[2], h = ee[3], k = ee[4], L = ee[5];
	double h_ = sqrt(p / mu);
	double n = h_ / (1.0 + f * cos(L) + g * sin(L));
	double s2 = 1.0 + h * h + k * k;
	double hh_kk = h * h - k * k;
	double hk2 = 2.0 * h * k;

	double r = n * h_ * mu;
	rv[0] = r / s2 * (cos(L) + hh_kk * cos(L) + hk2 * sin(L));
	rv[1] = r / s2 * (sin(L) - hh_kk * sin(L) + hk2 * cos(L));
	rv[2] = 2.0 * r / s2 * (h * sin(L) - k * cos(L));
	rv[3] = -1.0 / h_ / s2 * (sin(L) + hh_kk * sin(L) - hk2 * cos(L) + g - hk2 * f + hh_kk * g);
	rv[4] = -1.0 / h_ / s2 * (-cos(L) + hh_kk * cos(L) + hk2 * sin(L) - f + hk2 * g + hh_kk * f);
	rv[5] = 2.0 / h_ / s2 * (h * cos(L) + k * sin(L) + f * h + g * k);
	flag = 1;
	return;
}
void rv2ee(int& flag, double* ee, const double* RV, double mu = 3.98600441800e+14)
{
	flag = 0;
	if (mu <= 0.0)
		return;
	int i;
	double R[3] = { RV[0], RV[1], RV[2] };
	double V[3] = { RV[3], RV[4], RV[5] };
	double radius = V_Norm2(R, 3);//����
	double velocity = V_Norm2(V, 3);//�ٶ�	
	double unitR[3];
	for (i = 0; i < 3; i++) unitR[i] = R[i] / radius;//����λʸ��    
	double unitV[3];
	for (i = 0; i < 3; i++) unitV[i] = V[i] / velocity;//����λʸ��
	double hvector[3];
	V_Cross(hvector, unitR, unitV);
	double h = radius * velocity * V_Norm2(hvector, 3);//�Ƕ���ֵ
	double unith[3];
	for (i = 0; i < 3; i++) unith[i] = hvector[i] / V_Norm2(hvector, 3);//����淨��λʸ��

	// unith=[sin(i)*sin(OMEGA),
	//       -sin(i)*cos(OMEGA),
	//       cos(i)];
	//ƫ����ʸ��	
	double evector[3];
	V_Cross(evector, unitV, unith);
	for (i = 0; i < 3; i++) evector[i] = (velocity * h / mu) * evector[i] - unitR[i];
	//���ܾ�����ʸ��,ģΪe
	double qvector[3];
	V_Cross(qvector, unith, evector);
	//�����һ��ʸ��
	double unitA[3];
	for (i = 0; i < 3; i++) unitA[i] = h / mu * V[i] - qvector[i];
	//������ϵʽ
	// evector=e*[cos(omega)*cos(OMEGA)-sin(omega)*sin(OMEGA)*cos(i),
	//            cos(omega)*sin(OMEGA)+sin(omega)*cos(OMEGA)*cos(i),
	//            sin(omega)*sin(i)];
	// qvector=e*[-sin(omega)*cos(OMEGA)-cos(omega)*sin(OMEGA)*cos(i),
	//            -sin(omega)*sin(OMEGA)+cos(omega)*cos(OMEGA)*cos(i),
	//            cos(omega)*sin(i)];
	// unitR=[cos(OMEGA)*cos(omega+f)-cos(i)*sin(OMEGA)*sin(omega+f),
	//        sin(OMEGA)*cos(omega+f)+cos(i)*cos(OMEGA)*sin(omega+f),
	//        sin(i)*sin(omega+f)];
	// unitA=[-cos(OMEGA)*sin(omega+f)-cos(i)*sin(OMEGA)*cos(omega+f),
	//        -sin(OMEGA)*sin(omega+f)+cos(i)*cos(OMEGA)*cos(omega+f),
	//         sin(i)*cos(omega+f)];
	//�Ƶ��ó�
	//evector[0]+qvector[1]=(1+cos(i))*e*cos(omega+OMEGA);
	//evector[1]-qvector[0]=(1+cos(i))*e*sin(omega+OMEGA);
	//unith[0]=(1+cos(i))*tan(i/2)*sin(OMEGA);
	//unith[1]=-(1+cos(i))*tan(i/2)*cos(OMEGA);
	// unitR[0]+unitA[1]=(1+cos(i))*cos(omega+OMEGA+f);
	// unitR[1]-unitA[0]=(1+cos(i))*sin(omega+OMEGA+f);

	ee[0] = h * h / mu;//p

	if (unith[2] + 1.0 <= 0.0)
	{
		//		cout<<"�����ǽӽ�180�ȣ����ʺ�����������������. "<<endl;
		return;
	}
	double cosiadd1 = 1.0 + unith[2];
	ee[1] = (evector[0] + qvector[1]) / cosiadd1;//f
	ee[2] = (evector[1] - qvector[0]) / cosiadd1;//g
	ee[3] = -unith[1] / cosiadd1;//h
	ee[4] = unith[0] / cosiadd1;//k
	ee[5] = atan2(unitR[1] - unitA[0], unitR[0] + unitA[1]);//L
	//ת����[0,2pi)��
	ee[5] = fmod(ee[5], D2PI);
	if (ee[5] < 0.0)
		ee[5] += D2PI;
	flag = 1;
	return;
}
