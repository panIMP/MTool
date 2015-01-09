#include "MahaDist.h"
#include "../SharedInc/EquationSolve.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int calcMatCovar(const double*** matXs, double** matCovar, int dim, int num)
{
	int d = 0;
	int d1 = 0;
	int d2 = 0;
	int n = 0;
	double** matXEven = NULL;
	double** matDif = NULL;
	double** matDifRevert = NULL;
	if (callocMat(&matXEven, dim, 1) < 0)
		return -1;
	if (callocMat(&matDif, num, dim) < 0)
		return -1;
	if (callocMat(&matDifRevert, dim, num) < 0)
		return -1;

	for (n = 0; n < num; ++n)
	{
		for (d = 0; d < dim; ++d)
		{
			matXEven[d][0] += matXs[n][d][0];
		}
	}
	for (d = 0; d < dim; ++d)
	{
		matXEven[d][0] /= (double)num;
	}

	for (n = 0; n < num; ++n)
	{
		for (d = 0; d < dim; ++d)
		{
			matDif[n][d] = matXs[n][d][0] - matXEven[d][0];
		}
	}

	revertMat(matDif, matDifRevert, num, dim);

	mulMat(matDifRevert, matDif, matCovar, dim, n, dim);

	for (d1 = 0; d1 < dim; ++d1)
	{
		for (d2 = 0; d2 < dim; ++d2)
		{
			matCovar[d1][d2] /= (double)(num - 1);
		}
	}

	free(matXEven);
	matXEven = NULL;
	free(matDif);
	matDif = NULL;
	free(matDifRevert);
	matDifRevert = NULL;

	return 0;
}


int getMahaDistance(const double** mat1, const double** mat2, const double** matCovarInv, int dim, double* dist)
{
	int d = 0;
	double** matDif = NULL;
	double** matDifRevert = NULL;
	double** t_mat = NULL;
	if (callocMat(&matDif, dim, 1) < 0)
		return -1;
	if (callocMat(&matDifRevert, 1, dim) < 0)
		return -1;
	if (callocMat(&t_mat, 1, dim) < 0)
		return -1;

	for (d = 0; d < dim; ++d)
	{
		matDif[d][0] = mat1[d][0] - mat2[d][0];
	}
	revertMat(matDif, matDifRevert, dim, 1);

	mulMat(matDifRevert, matCovarInv, t_mat, 1, dim, dim);
	mulMat(t_mat, matDif, &dist, 1, dim, 1);
	*dist = sqrt(*dist);

	free(matDif);
	matDif = NULL;
	free(matDifRevert);
	matDifRevert = NULL;
	free(t_mat);
	t_mat = NULL;

	return 0;
}


int getMahaDistance2(const double** mat1, const double** mat2, const double** matCovarInv, int dim, double* dist)
{
	int d = 0;
	double** matDif = NULL;
	double** matDifRevert = NULL;
	double** t_mat = NULL;
	if (callocMat(&matDif, dim, 1) < 0)
		return -1;
	if (callocMat(&matDifRevert, 1, dim) < 0)
		return -1;
	if (callocMat(&t_mat, 1, dim) < 0)
		return -1;

	for (d = 0; d < dim; ++d)
	{
		matDif[d][0] = mat1[d][0] - mat2[d][0];
	}
	revertMat(matDif, matDifRevert, dim, 1);

	mulMat(matDifRevert, matCovarInv, t_mat, 1, dim, dim);
	mulMat(t_mat, matDif, &dist, 1, dim, 1);

	free(matDif);
	matDif = NULL;
	free(matDifRevert);
	matDifRevert = NULL;
	free(t_mat);
	t_mat = NULL;

	return 0;
}


/*
// sample
int main()
{
	int num = 0;
	int dim = 0;
	int n = 0;
	int d = 0;
	int d2 = 0;
	double*** matXs = NULL;
	double** matCovar = NULL;
	double** matCovarInv = NULL;
	double dist = 0.0;

	fprintf(stdout, "input the dim of one vector: ");
	fscanf_s(stdin, "%d", &dim);
	fprintf(stdout, "input the number of vectors: ");
	fscanf_s(stdin, "%d", &num);

	fprintf(stdout, "input the value of each element: \n");
	matXs = (double***)malloc(num * sizeof(double**));
	for (n = 0; n < num; ++n)
	{
		if (mallocMat(&matXs[n], dim, 1) < 0)
			return -1;

		for (d = 0; d < dim; ++d)
		{
			fscanf_s(stdin, "%lf", &matXs[n][d][0]);
		}
	}

	if (mallocMat(&matCovar, dim, dim) < 0)
		return -1;
	
	calcMatCovar(matXs, matCovar, dim, num);

	for (d = 0; d < dim; ++d)
	{
		for (d2 = 0; d2 < dim; ++d2)
		{
			fprintf(stdout, "%lf\t", matCovar[d][d2]);
		}
		fprintf(stdout, "\n");
	}

	if (mallocMat(&matCovarInv, dim, dim) < 0)
		return -1;
	if (inverseMat(matCovar, matCovarInv, dim) < 0)
		return -1;

	for (d = 0; d < dim; ++d)
	{
		for (d2 = 0; d2 < dim; ++d2)
		{
			fprintf(stdout, "%lf\t", matCovarInv[d][d2]);
		}
		fprintf(stdout, "\n");
	}

	if (getMahaDistance(matXs[1], matXs[3], matCovarInv, dim, &dist) < 0)
		return -1;
	fprintf(stdout, "maha dist of vector 0 and 1 is:%lf", dist);

	free(matXs);
	matXs = NULL;
	free(matCovar);
	matCovar = NULL;
	free(matCovarInv);
	matCovarInv = NULL;

	system("pause");

	return 0;
}
*/