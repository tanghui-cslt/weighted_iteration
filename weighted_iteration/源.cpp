#include <iostream>
#include <iomanip>
using namespace std;
#define MAX_TIMES 100
#define ERROR 1e-11

void scanf_data(double ** &A, double *&x, double *&b, int &n);
void LU_decom(double **&a, double **&l, double **&u, int n);
double * Calc_X(double **&l, double **&u, double *&b, int n);
void weighted_iteration(double ** &A, double *&x, double *&b, int &n);
double * Calc_r(double ** &A, double *&x, double *&b, int &n);
double *sum_vector(double *a, double *b, int n);
double max_element(double *a, int n);

int main()
{
	double **A = NULL, *x = NULL, *b = NULL;
	int n = 0;

	scanf_data(A, x, b, n);

	weighted_iteration(A, x, b, n);

	getchar();
	getchar();
	return 0;
}
double calc_alpha(double **A, const int n)
{
	double alpha, min = A[0][0];
	for (int i = 0; i< n; i ++)
	{
		for (int j = 0; j < n; j++)
		{
			min = min > A[i][j]? A[i][j] : min;
		}
	}
	//cout << "min = " << min;
	alpha = min*min / 9;
	return alpha;
}
void Correction_A(double **&A, double alpha, const int n)
{
	for (int i = 0; i< n; i++)
	{
		A[i][i] += alpha;
	}
}
void weighted_iteration(double ** &A, double *&x, double *&b, int &n)
{
	double **l = NULL, **u = NULL;
	double *r = NULL, *d = NULL;
	double norm_x = 0, norm_d = 0;
	int i = 0;
	double alpha = calc_alpha(A, n);
	//alpha = 0.00011;
	//cout << "alpha = " << alpha << endl;
	
	Correction_A(A,alpha,n);


	LU_decom(A, l, u, n);


	x = Calc_X(l, u, b, n);
	
	for (i = 0; i < MAX_TIMES; i++)
	{

		r = Calc_r(A, x, b, n);

		d = Calc_X(l, u, r, n);

		x = sum_vector(x, d, n);

		norm_x = max_element(x, n);

		norm_d = max_element(d, n);

		double min = fabs(norm_d / norm_x);
	//	cout << "min = " << min << endl;
		if (min < ERROR)
		{
			cout << "达到精确值\n";
			break;
		}
	}
	cout << "迭代次数：" << i << endl;
	if (i == MAX_TIMES)
	{
		cout << "达到最大迭代次数\n";
	}
	cout << "x的值为：\n";
	for (int i = 0; i < n; i++)
	{
		cout << "\tx" << "[" << i << "]=" << x[i];
	}
	cout << endl;
}

double max_element(double *a, int n)
{
	double max_num = a[0];
	for (int i = 1; i < n; i++)
	{
		max_num = max_num < a[i] ? a[i] : max_num;
	}
	return max_num;
}
double *sum_vector(double *a, double *b, int n)
{
	for (int i = 0; i < n; i++)
	{
		a[i] += b[i];
	}
	return a;
}

double * Calc_r(double ** &A, double *&x, double *&b, int &n)
{
	double *r = new double[n];

	for (int i = 0; i < n; i++)
	{
		double sum = 0;
		for (int j = 0; j < n; j++)
		{
			sum += A[i][j] * x[j];
		}
		r[i] = b[i] - sum;
		//cout << r[i] << " ";
	}
	return r;
}
void scanf_data(double ** &A, double *&x, double *&b, int &n)
{
	cout << "请输入矩阵的维数:\n";
	cin >> n;

	A = new double *[n];
	b = new double[n];
	cout << "请输入系数矩阵A:\n";
	for (size_t i = 0; i < n; i++)
	{
		A[i] = new double[n];
		for (size_t j = 0; j < n; j++)
		{
			A[i][j] = 1.0/(i+1+j+1);
		}
	}

	//x = new double[n];
	//cout << "请输入x初值:\n";
	//for (size_t i = 0; i < n; i++)
	//{
	//	cin >> x[i];
	//}

	cout << "请输入向量b:\n";
	for (size_t i = 0; i < n; i++)
	{
		b[i] = 0;
		for (int j = 0; j < n; j++)
		{
			b[i] += A[i][j];
		}
		//cin >> b[i];
		
	}
}

void LU_decom(double **&a, double **&l, double **&u, int n)
{
	//初始化空间l u
	l = new double *[n];
	u = new double *[n];
	for (int i = 0; i < n; i++)
	{
		l[i] = new double[n];
		u[i] = new double[n];
		for (int j = 0; j < n; j++)
		{
			l[i][j] = 0;
			u[i][j] = 0;
			if (i == j)
				l[i][j] = 1;
		}
	}

	//计算l u矩阵的结果
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (j == 0 && i != 0)
				l[i][j] = a[i][j] / u[0][0];
			else if (i == 0)
				u[i][j] = a[i][j];
			else {
				int min_temp = i > j ? j : i;

				if (i > j)
				{
					double sum = 0;
					for (int k = 0; k < min_temp; k++)
					{
						sum += l[i][k] * u[k][j];
					}
					l[i][j] = (a[i][j] - sum) / u[min_temp][min_temp];
				}
				else {
					double sum = 0;
					for (int k = 0; k < min_temp; k++)
					{
						sum += l[i][k] * u[k][j];
					}
					u[i][j] = (a[i][j] - sum);
				}
			}
			//cout << l[i][j] << " ";
		}
		//cout << endl;
	}

}

double * Calc_X(double **&l, double **&u, double *&b, int n)
{
	double *y = NULL, *x = NULL;
	x = new double[n];
	y = new double[n];
	//计算y
	for (int i = 0; i < n; i++)
	{

		double sum_x = 0;
		for (int j = 0; j < i; j++)
		{
			sum_x += l[i][j] * y[j];
		}
		y[i] = b[i] - sum_x;
		//cout << y[i] << " ";
	}

	//计算x

	for (int i = n - 1; i >= 0; i--)
	{

		double sum_y = 0;
		for (int j = n - 1; j > i; j--)
		{
			sum_y += u[i][j] * x[j];
		}
		x[i] = (y[i] - sum_y) / u[i][i];
	}


	return x;
}