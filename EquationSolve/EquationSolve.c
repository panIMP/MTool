#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_MATRIX_LEN 10


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
		det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];

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



int main()
{
	int n;
	int r, c;
	double** matA = NULL;
	double** matAInv = NULL;
	double** matB = NULL;
	double** matX = NULL;

	printf("input the n of matrix: ");
	fscanf_s(stdin, "%d", &n);
	
	mallocMat(&matA, n, n);
	mallocMat(&matB, n, 1);
	mallocMat(&matAInv, n, n);

	printf("input the value of the matA: \n");
	for (r = 0; r < n; ++r)
	{
		for (c = 0; c < n; ++c)
		{
			fscanf_s(stdin, "%lf", &matA[r][c]);
		}
	}

	inverseMat(matA, matAInv, n);
	for (r = 0; r < n; ++r)
	{
		for (c = 0; c < n; ++c)
		{
			printf("%lf\t", matAInv[r][c]);
		}
		printf("\n");
	}

	getchar();

	return 0; 
}