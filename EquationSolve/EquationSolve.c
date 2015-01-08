#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_MATRIX_LEN 10
#define EPS 1e-30


int mallocMat(double*** matA, int rows, int cols)
{
	int i, j;

	*matA = (double**)malloc(sizeof(double*) * rows);
	if (*matA == NULL)
		return -1;

	for (i = 0; i < rows; ++i)
	{
		if (((*matA)[i] = (double*)malloc(sizeof(double)* cols)) == NULL)
		{
			return -1;
		}
	}

	return 0;
}


int calcMatDet(double** mat, double* det, int n)
{
	int r = 0; 
	int c = 0;
	double** t_mat = NULL;
	int t_r = 0;
	int t_c = 0;

	if (n == 2)
	{
		*det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
		return 0;
	}

	if (mallocMat(&t_mat, n - 1, n - 1) < 0)
	{
		return -1;
	}

	for (r = 0; r < n; ++r)
	{
		int w = pow(-1, r);
		double t_det = 0.0;

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

		if (calcMatDet(t_mat, &t_det, n - 1) < 0)
		{
			return -1;
		}
		(*det) += w * mat[r][0] * t_det;
	}

	free(t_mat);
	t_mat = NULL;

	return 0;
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
	double** matRev = NULL;
	int r = 0;
	int c = 0;
	int t_r = 0;
	int t_c = 0;
	double det = 0.0;

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
		return -1;
	}
	if (mallocMat(&matRev, n, n) < 0)
	{
		return -1;
	}

	for (r = 0; r < n; ++r)
	{
		for (c = 0; c < n; ++c)
		{
			int w = pow(-1, (r + c) % 2);
			double t_det = 0.0;

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

			if (calcMatDet(t_mat, &t_det, n - 1) < 0)
			{
				return -1;
			}
			matRev[r][c] = w * t_det;
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

	free(matRev);
	matRev = NULL;
	free(t_mat);
	t_mat = NULL;

	return 0;
}


int mulMat(double** mat1, double** mat2, double** matRes, int rows1, int cols1rows2, int cols2)
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
			for (c1 = 0; c1 < cols1rows2; ++c1)
			{
				r2 = c1;
				val += mat1[r1][c1] * mat2[r2][c2];
			}

			matRes[r1][c2] = val;
			val = 0.0;
		}
	}

	return 0;
}


int calcEquationSolution(double** matA, double** matB, double** matX, int n)
{
	double** matAInv = NULL;

	if (mallocMat(&matAInv, n, n) < 0)
		return -1;

	if (inverseMat(matA, matAInv, n) < 0)
	{
		return -1;
	}

	mulMat(matAInv, matB, matX, n, n, 1);

	free(matAInv);
	matAInv = NULL;

	return 0;
}


int calcMatTransformation(double** src, double** dst, double** matT, int n)
{
	double** srcInv = NULL;

	if (mallocMat(&srcInv, n, n) < 0)
		return -1;

	if (inverseMat(src, srcInv, n) < 0)
		return -1;
	
	mulMat(dst, srcInv, matT, n, n, n);

	free(srcInv);
	srcInv = NULL;

	return 0;
}


int main()
{
	int n;
	int r, c;
	double** matA = NULL;
	double** matX = NULL;
	double** matB = NULL;

	fprintf(stdout, "input the n of matrix: ");
	fscanf_s(stdin, "%d", &n);
	mallocMat(&matA, n, n);
	mallocMat(&matB, n, n);
	mallocMat(&matX, n, n);

	fprintf(stdout, "input the value of the matA: \n");
	for (r = 0; r < n; ++r)
	{
		for (c = 0; c < n; ++c)
		{
			fscanf_s(stdin, "%lf", &matA[r][c]);
		}
	}

	fprintf(stdout, "input the value of the matB: \n");
	for (r = 0; r < n; ++r)
	{
		for (c = 0; c < n; ++c)
		{
			fscanf_s(stdin, "%lf", &matB[r][c]);
		}
	}

	calcMatTransformation(matA, matB, matX, n);
	for (r = 0; r < n; ++r)
	{
		for (c = 0; c < n; ++c)
		{
			fprintf(stdout, "%lf\t", matX[r][c]);
		}
	}

	system("pause");
	return 0; 
}