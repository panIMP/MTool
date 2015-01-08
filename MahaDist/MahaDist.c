#include "MahaDist.h"
#include "../SharedInc/EquationSolve.h"

#include <stdio.h>
#include <stdlib.h>


int calcMatCovar(double*** matXs, double** matCovar, int dim, int num)
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



int main()
{
	int num = 0;
	int dim = 0;
	int n = 0;
	int d = 0;
	int d2 = 0;
	double*** matXs = NULL;
	double** matCovar = NULL;

	fprintf(stdout, "input the value of vector dim: ");
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

	if (inverseMat(matCovar, matCovar, dim) < 0)
		return -1;

	for (d = 0; d < dim; ++d)
	{
		for (d2 = 0; d2 < dim; ++d2)
		{
			fprintf(stdout, "%lf\t", matCovar[d][d2]);
		}
		fprintf(stdout, "\n");
	}

	free(matXs);
	matXs = NULL;

	system("pause");

	return 0;
}