#include <iostream>
#include <iomanip>
#include <time.h>
#include <io.h>
#include <omp.h>

double ** Generate_Matrix(int N)
{
	double ** Matrix;
	Matrix = new double *[N];
	for (int i = 0; i < N; i++)
		Matrix[i] = new double[N];
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
		Matrix[i][j] = 0 + (100 - 0)* ((double)rand() / RAND_MAX);
	return Matrix;
}

void Basic_Matrix_Multiplication(double **MatrixC, double **MatrixA, double **MatrixB, int N)
{
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
	for (int k = 0; k < N; k++)
		MatrixC[i][j] = MatrixC[i][j] + MatrixA[i][k] * MatrixB[k][j];
}

void Matrix_Multiplication_Change_jk(double **MatrixC, double **MatrixA, double **MatrixB, int N)
{
	for (int i = 0; i < N; i++)
	for (int k = 0; k < N; k++)
	for (int j = 0; j < N; j++)
		MatrixC[i][j] += MatrixA[i][k] * MatrixB[k][j];
}


void MatrixMultiplication_Vectorization(double **MatrixC, double **MatrixA, double **MatrixB, int N)
{
	for (int i = 0; i < N; i++)
	for (int k = 0; k < N; k++)
#pragma ivdep
	for (int j = 0; j < N; j++)
		MatrixC[i][j] = MatrixC[i][j] + MatrixA[i][k] * MatrixB[k][j];
}


void MatrixMultiplication_omp_parallel(double **MatrixC, double **MatrixA, double **MatrixB, int size)
{
#pragma omp parallel for
	for (int i = 0; i < size; i++)
	{
		for (int k = 0; k < size; k++)
#pragma ivdep
		for (int j = 0; j < size; j++)
		{
			MatrixC[i][j] = MatrixC[i][j] + MatrixA[i][k] * MatrixB[k][j];
		}
	}
}

void ZeroC(double **MatrixC, int N)//обнуление матрицы
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			MatrixC[i][j] = 0;
	}
}

void Matrix_Print(double **MatrixA, int N)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			std::cout << "  " << MatrixA[i][j];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

int Validation_Check(double** Matrix1, double** Matrix2, int N)
{
	double eps = 0.001;
	double Difference_max = abs(Matrix1[0][0] - Matrix2[0][0]);
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
	{
		if (abs(Matrix1[i][j] - Matrix2[i][j]) > Difference_max)
			Difference_max = abs(Matrix1[i][j] - Matrix2[i][j]);
	}
	if (Difference_max>eps)
		return 0;
	else return 1;
}

int main(int argc, char* argv[])
{
	int N = 2200;
	int number = 0;
	double **MatrixA;
	double **MatrixB;
	double **MatrixC;
	double **MatrixC_basic;
	MatrixA = Generate_Matrix(N);
	MatrixB = Generate_Matrix(N);
	MatrixC = new double *[N];
	for (int i = 0; i < N; i++)
	{
		MatrixC[i] = new double[N];
	}
	MatrixC_basic = new double *[N];
	for (int i = 0; i < N; i++)
	{
		MatrixC_basic[i] = new double[N];
	}
	clock_t start, end_;
	if (argc != 2)
		std::cout << "Enter argument!" << std::endl;
	else
		number = atoi(argv[1]);
	анулируем время 
	srand((unsigned)time(NULL));
	//ZeroC(MatrixC_basic, N);
	ZeroC(MatrixC, N);
	if (number == 1)
	{
		start = clock();
		//Basic_Matrix_Multiplication(MatrixC_basic, MatrixA, MatrixB, N);
		Basic_Matrix_Multiplication(MatrixC, MatrixA, MatrixB, N);
		end_ = clock();
		printf("Function Basic_Matrix_Multiplication was executed in %.4f second(s)\n", ((double)end_ - start) / ((double)CLOCKS_PER_SEC));
	}
	else if (number == 2)
		{
			start = clock();
			Matrix_Multiplication_Change_jk(MatrixC, MatrixA, MatrixB, N);
			end_ = clock();
			printf("Function Matrix_Multiplication_Change_jk was executed in %.4f second(s)\n", ((double)end_ - start) / ((double)CLOCKS_PER_SEC));
			/*if (Validation_Check(MatrixC_basic, MatrixC, N)==0)
				printf("Not correct!!!\n");*/
		}
	else if (number == 3)
	{
		start = clock();
		MatrixMultiplication_Vectorization(MatrixC, MatrixA, MatrixB, N);
		end_ = clock();
		printf("Function MatrixMultiplication_Vectorization was executed in %.4f second(s)\n", ((double)end_ - start) / ((double)CLOCKS_PER_SEC));
		/*if (Validation_Check(MatrixC_basic, MatrixC, N) == 0)
			printf("Not correct!!!\n");*/
	}

	else if (number == 4)
	{
		start = clock();
	MatrixMultiplication_omp_parallel(MatrixC, MatrixA, MatrixB, N);
	end_ = clock();
	printf("Function MatrixMultiplication_omp_parallel was executed in %.4f second(s)\n", ((double)end_ - start) / ((double)CLOCKS_PER_SEC));
	/*if (Validation_Check(MatrixC_basic, MatrixC, N) == 0)
		printf("Not correct!!!\n");*/
	}
	for (int i = 0; i < N; i++)
	{
		delete[] MatrixA[i];
		delete[] MatrixB[i];
		delete[] MatrixC[i];
		delete[] MatrixC_basic[i];
	}
	delete[] MatrixA;
	delete[] MatrixB;
	delete[] MatrixC;
	delete[] MatrixC_basic;
}



