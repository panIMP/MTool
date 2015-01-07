#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_MATRIX_LEN 10
#define EPS 1e-6


int mallocMat(double*** matA, int rows, int cols)
{
	int i, j;

	*matA = (double**)malloc(sizeof(double*) * rows);
	if (*matA == NULL)
		return -1;

	for (i = 0; i < rows; ++i)
	{
		if (((*matA)[i] = (double*)malloc(sizeof(double)* cols)) == NULL)
			return -1;
	}

	return 0;
}


double calcDetMat(double** mat, int n)
{
	double det = 0.0;
	int r = 0; 
	int c = 0;
	double** t_mat = NULL;
	int t_r = 0;
	int t_c = 0;

	if (n == 2)
	{
		det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
		return det;
	}

	if (mallocMat(&t_mat, n - 1, n - 1) < 0)
	{
		printf("malloc failure in mallocMat!");
		exit(-1);
	}

	for (r = 0; r < n; ++r)
	{
		int w = pow(-1, r);

		for (t_r = 0; t_r < n - 1; ++t_r)
		{
			for (t_c = 0; t_c < n - 1; ++t_c)
			{
				int rCur = t_r;

				if (t_r >= r)
					rCur++;

				t_mat[t_r][t_c] = mat[rCur][t_c + 1];
			}
		}

		det += w * mat[r][0] * calcDetMat(t_mat, n - 1);
	}

	return det;
}


int revertMat(double** mat, double** matRev, int rows, int cols)
{
	int r = 0; 
	int c = 0;

	for (r = 0; r < rows; ++r)
	{
		for (c = 0; c < cols; ++c)
		{
			matRev[c][r] = mat[r][c];
		}
	}

	return 0;
}


int inverseMat(double** mat, double** matInv, int n)
{
	double** t_mat = NULL;
	int r = 0;
	int c = 0;
	int t_r = 0;
	int t_c = 0;
	double det = 0.0;
	double** matRev = NULL;

	if (n == 1)
	{
		matInv[0][0] = 1.0 / mat[1][1];
		
		return 0;
	}

	if (n == 2)
	{
		if ((det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) == 0.0)
			return -1;

		matInv[0][0] =  mat[1][1] / det;
		matInv[0][1] = -mat[0][1] / det;
		matInv[1][0] = -mat[1][0] / det;
		matInv[1][1] =  mat[0][0] / det;

		return 0;
	}

	if (mallocMat(&t_mat, n-1, n-1) < 0)
	{
		printf("malloc failure in mallocMat!");
		return -1;
	}
	if (mallocMat(&matRev, n, n) < 0)
	{
		printf("malloc failure in mallocMat!");
		return -1;
	}

	for (r = 0; r < n; ++r)
	{
		for (c = 0; c < n; ++c)
		{
			int w = pow(-1, (r + c) % 2);

			for (t_r = 0; t_r < n - 1; ++t_r)
			{
				for (t_c = 0; t_c < n - 1; ++t_c)
				{
					int rCur = t_r;
					int cCur = t_c;

					if (t_r >= r)
						rCur++;

					if (t_c >= c)
						cCur++;

					t_mat[t_r][t_c] = mat[rCur][cCur];
				}
			}

			matRev[r][c] = w * calcDetMat(t_mat, n - 1);
		}

		det += matRev[r][0] * mat[r][0];
	}

	if (det == 0.0)
		return -1;

	for (r = 0; r < n; ++r)
	{
		for (c = 0; c < n; ++c)
		{
			matRev[r][c] /= det;
		}
	}

	revertMat(matRev, matInv, n, n);

	return 0;
}


int mulMat(double** mat1, double** mat2, double** mat3, int rows1, int n, int cols2)
{
	int r1 = 0;
	int c1 = 0;
	int r2 = 0;
	int c2 = 0;
	double val = 0.0;

	for (r1 = 0; r1 < rows1; ++r1)
	{
		for (c2 = 0; c2 < cols2; ++c2)
		{
			for (c1 = 0; c1 < n; ++c1)
			{
				r2 = c1;
				val += mat1[r1][c1] * mat2[r2][c2];
			}

			mat3[r1][c2] = val;
			val = 0.0;
		}
	}

	return 0;
}


int main()
{
	int n;
	int r, c;
	double** matA = NULL;
	double** matAInv = NULL;
	double** matB = NULL;
	double** matX = NULL;

	// input matA, matB
	printf("input the n of matrix: ");
	fscanf_s(stdin, "%d", &n);
	mallocMat(&matA, n, n);
	mallocMat(&matB, n, 1);
	mallocMat(&matAInv, n, n);
	mallocMat(&matX, n, n);

	printf("input the value of the matA: \n");
	for (r = 0; r < n; ++r)
	{
		for (c = 0; c < n; ++c)
		{
			fscanf_s(stdin, "%lf", &matA[r][c]);
		}
	}

	printf("input the value of matB: \n");
	for (r = 0; r < n; ++r)
	{
		fscanf_s(stdin, "%lf", &matB[r][0]);
	}

	// get inverse mat of matA and output inverse mat of matA
	if (inverseMat(matA, matAInv, n) < 0)
	{
		printf("inverse mat not exists");
		return -1;
	}

	printf("the inverse mat of matA is: \n");
	for (r = 0; r < n; ++r)
	{
		for (c = 0; c < n; ++c)
		{
			printf("%lf\t", matAInv[r][c]);
		}
		printf("\n");
	}

	// get matX and output matX
	mulMat(matAInv, matB, matX, n, n, 1);
	printf("the solution is: \n");
	for (r = 0; r < n; ++r)
	{
		printf("%lf\t", matX[r][0]);
	}

	getchar();
	getchar();
	return 0; 
}