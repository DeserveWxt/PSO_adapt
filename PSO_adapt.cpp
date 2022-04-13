#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<iostream>
#include<algorithm>
#include<iomanip>

using namespace std;

#define CFYZ	0.1				// 惩罚因子			0.1， 1， 10， 100……
#define RWSJ	100				// 任务数据大小		100 - 150 Kb
#define RWQD	2000				// 任务强度			2000-3000 cycles/bit

#define M	20				// 粒子数目
#define K	10				// 用户设备数目
#define D	4				// 解空间维数
#define Pmax	0.6				// 用户设备的最大传输功率		600 mw
#define F	(8*pow(10,10))			// MEC服务器的最大计算能力	8 * 10^10 cycles/s
#define Wmax	0.8				// 最大惯性权重
#define Wmin	0.4				// 最小惯性权重
#define C1	2				// 学习因子
#define C2	2				// 学习因子
#define e	(pow(10,-6))			// 收敛因子
#define T	1000				// 最大迭代周期
#define B	(4*pow(2,20))			// 信道带宽			4 MHz
#define N0	(pow(10,-14))			// 噪声功率谱密度		10^-11 mV/Hz
#define SJ	(((double)rand()) / RAND_MAX)	// 0-1之间的随机数		// (((double)(rand() % 99) + 1) / 100)

double TotalFitness;				// 每次迭代总适应度值
double FavFitness;				// 每次迭代平均适应度值
double MaxFitness;				// 每次迭代最优适应度值
double MinExeTime;				// 总执行时延
double TotalPenalty;

double P[M][K][D];				// 种群位置数组
double V[M][K][D];				// 种群速度数组
double Pbest[M][K][D];				// 种群粒子个体极值
double Gbest[K][D];				// 种群粒子全局极值

double E[K];					// 惩罚因子数组
double Data[K];					// 任务大小数组				100  -  150  kb
double C[K];					// 任务强度，计算1bit所需CPU周期数		2000  -  3000  cycles/bit

double H[K] =      { 1.3,  1.6,  1.4,  1.8,  1.2,  1.7,  1.5,  1.1,  1.9,   1  };	// 用户设备信道增益数组
double J[K] =      { 200,  240,  280,  290,  230,  210,  250,  270,  220,  260 };	// 压缩1bit所需CPU周期数				200  -  500  cycles/bit
double L[K] =      { 200,  240,  280,  290,  230,  210,  250,  270,  220,  260 };	// 解压1bit所需CPU周期数				200  -  500  cycles/bit
double Energy[K] = { 1.7,  1.9,  1.5,  1.6,  1.5,   2,   1.9,  1.8,  1.7,  1.6 };	// 用户设备CPU每运行一周期的能耗		1.5*10^-11  -  2.0*10^-11  J/cycle
double b[K] =      { 0.5,  0.4,  0.7,  0.5,  0.6,  0.4,  0.7,  0.6,  0.5,  0.4 };	// 压缩系数数组							{ 0.4, 0.5, 0.6, 0.7 }
double LE[K] =     { 1.3,  1.2,  1.5,   1,   1.4,  1.2,  1.5,  1.3,  1.5,   1  };	// 用户设备本地计算频率数组				{ 1, 1.1, 1.2, 1.3, 1.4, 1.5 }GHz

double ExeTime[K];				// 任务执行时间
double ExeEnergy[K];				// 任务执行能耗
double Penalty[K];				// 

double TotalTime[M];				//
double W[M];					// 种群粒子惯性权重数组
double Fitness[M];				// 种群粒子适应度数组
double FitnessBest[M];				// 种群粒子最优适应度数组
double FitnessBest2[2];				// 种群两代最优适应度

double dmax(double a, double b)
{
	if (a > b)
		return a;
	else
		return b;
}

int main()
{
	srand(time(0));

	// 初始化
	for (int i = 0; i < K; i++)
	{
		E[i] = CFYZ;				// 惩罚因子
		Data[i] = RWSJ * pow(2, 10);		// 任务数据大小
		C[i] = RWQD;				// 任务强度
		Energy[i] = Energy[i] * pow(10, -11);	// 用户设备CPU每运行一周期的能耗
		LE[i] = LE[i] * pow(2, 30);		// 用户设备本地计算频率
	}

	//初始化种群粒子位置
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < K; j++)
		{
			P[i][j][0] = SJ;		// 任务卸载比例
			P[i][j][1] = SJ;		// 卸载任务中选取压缩部分的比例
			P[i][j][2] = SJ * Pmax;		// 用户设备的传输功率
			P[i][j][3] = SJ * F;		// MEC服务器分配给用户设备的计算资源
		}
	}

	// 对MEC分配给用户设备的计算资源进行归一化调整
	for (int i = 0; i < M; i++)
	{
		double TotalF = 0;
		for (int j = 0; j < K; j++)
		{
			TotalF += P[i][j][3];
		}
		for (int j = 0; j < K; j++)
		{
			P[i][j][3] = (P[i][j][3] / TotalF) * F;
		}
	}

	// 初始化种群粒子速度
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < K; j++)
		{
			V[i][j][0] = (-0.2 + 0.4 * SJ) * SJ;
			V[i][j][1] = (-0.2 + 0.4 * SJ) * SJ;
			V[i][j][2] = (-0.2 + 0.4 * SJ) * SJ * Pmax;
			V[i][j][3] = (-0.2 + 0.4 * SJ) * SJ * P[i][j][3];
		}
	}

	// 开始迭代
	for (int t = 0; t < T; t++)
	{
		TotalFitness = 0;

		// 粒子个体适应度值计算及自适应惯性权重的设置
		for (int i = 0; i < M; i++)
		{
			TotalTime[i] = 0;
			TotalPenalty = 0;

			for (int j = 0; j < K; j++)
			{
				double tjl0 = Data[j] * P[i][j][0] * P[i][j][1] * J[j] / LE[j];								// 本地压缩时延
				double tjl1 = Data[j] * (1 - P[i][j][0]) * C[j] / LE[j];								// 本地计算时延
				double tjt0 = Data[j] * P[i][j][0] * (1 - P[i][j][1]) / (B * (log(1 + P[i][j][2] * H[j] / (N0 * B)) / log(2)));		// 卸载任务中未选取压缩任务的传输时延
				double tjt1 = Data[j] * P[i][j][0] * P[i][j][1] * b[j] / (B * (log(1 + P[i][j][2] * H[j] / (N0 * B)) / log(2)));	// 卸载任务中选取压缩任务的传输时延
				double tjc0 = Data[j] * P[i][j][0] * (1 - P[i][j][1]) * C[j] / P[i][j][3];						// 卸载任务中未选取压缩任务的MEC服务器的执行时延
				double tjc1 = Data[j] * P[i][j][0] * P[i][j][1] * (C[j] + L[j]) / P[i][j][3];						// 卸载任务中选取压缩任务的MEC服务器的执行时延
				ExeTime[j] = dmax(tjl0 + tjl1, dmax(tjt0 + tjc0, dmax(tjl0, tjt0) + tjt1) + tjc1);					// 该任务执行时延

				TotalTime[i] += ExeTime[j];		// 总时延

				double ejl0 = Data[j] * P[i][j][0] * P[i][j][1] * J[j] * Energy[j];		// 本地压缩能耗
				double ejl1 = Data[j] * (1 - P[i][j][0]) * C[j] * Energy[j];			// 本地计算能耗
				double ejt0 = P[i][j][2] * tjt0;						// 卸载任务中未选取压缩任务的传输能耗
				double ejt1 = P[i][j][2] * tjt1;						// 卸载任务中选取压缩任务的传输能耗
				ExeEnergy[j] = ejl0 + ejl1 + ejt0 + ejt1;					// 该任务执行能耗

				Penalty[j] = E[j] * dmax(0, ExeEnergy[j] - Data[j] * C[j] * Energy[j] * 0.9);

				TotalPenalty += Penalty[j];
			}

			Fitness[i] = -(TotalTime[i] + TotalPenalty);

			TotalFitness += Fitness[i];
		}

		FavFitness = TotalFitness / M;

		MaxFitness = *max_element(Fitness, Fitness + M);

		MinExeTime = *min_element(TotalTime, TotalTime + M);

		FitnessBest2[t % 2] = MaxFitness;

		for (int i = 0; i < M; i++)
		{
			if (Fitness[i] > FavFitness)
				W[i] = Wmin + (Wmax - Wmin) * (MaxFitness - Fitness[i]) / (MaxFitness - FavFitness);
			else
				W[i] = Wmax;
		}

		// 种群粒子最优适应度值初始化
		if (t == 0)
		{
			for (int i = 0; i < M; i++)
			{
				FitnessBest[i] = Fitness[i];
			}
		}

		// 种群粒子最优适应度值及个体极值的更新
		for (int i = 0; i < M; i++)
		{
			if (Fitness[i] >= FitnessBest[i])
			{
				FitnessBest[i] = Fitness[i];
				for (int j = 0; j < K; j++)
				{
					for (int k = 0; k < D; k++)
					{
						Pbest[i][j][k] = P[i][j][k];
					}
				}
			}
		}

		// 种群最优适应度值及全局极值的更新
		for (int i = 0; i < M; i++)
		{
			if (Fitness[i] == MaxFitness)
			{
				for (int j = 0; j < K; j++)
				{
					for (int k = 0; k < D; k++)
					{
						Gbest[j][k] = P[i][j][k];
					}
				}
				break;
			}
		}

		// 更新粒子的速度和位置
		for (int i = 0; i < M; i++)
		{
			for (int j = 0; j < K; j++)
			{
				for (int k = 0; k < D; k++)
				{
					V[i][j][k] = W[i] * V[i][j][k] + C1 * SJ * (Pbest[i][j][k] - P[i][j][k]) + C2 * SJ * (Gbest[j][k] - P[i][j][k]);
					P[i][j][k] = P[i][j][k] + V[i][j][k];
				}
			}
		}

		// 限制
		for (int i = 0; i < M; i++)
		{
			for (int j = 0; j < K; j++)
			{
				// 任务卸载比例
				if (P[i][j][0] >= 1)	P[i][j][0] = 0.99;
				if (P[i][j][0] <= 0)	P[i][j][0] = 0.01;

				// 卸载任务中选取压缩部分的比例
				if (P[i][j][1] >= 1)	P[i][j][1] = 0.99;
				if (P[i][j][1] <= 0)	P[i][j][1] = 0.01;

				// 用户设备的传输功率
				if (P[i][j][2] >= Pmax)	P[i][j][2] = Pmax;
				if (P[i][j][2] <= 0)	P[i][j][2] = 0.1 * Pmax;

				// MEC服务器分配给用户设备的计算资源
				if (P[i][j][3] >= F)	P[i][j][3] = 0.99 * F;
				if (P[i][j][3] <= 0)	P[i][j][3] = 0.01 * F;
			}
		}

		cout << "t=" << t + 1 << "		";
		cout << "MaxFitness = " << setw(12) << std::left << MaxFitness << "		";
		cout << "MinExeTime = " << setw(10) << std::left << MinExeTime << endl;

		//if (t != 0 && fabs(FitnessBest2[0] - FitnessBest2[1]) < e)	break;
	}

	return 0;
}
