#ifndef _EQUATION_SOLVE_H
#define _EQUATION_SOLVE_H

int mallocMat(double*** matA, int rows, int cols);

int callocMat(double*** matA, int rows, int cols);

int calcMatDet(const double** mat, double* det, int n);

int revertMat(const double** mat, double** matRev, int rows, int cols);

int inverseMat(const double** mat, double** matInv, int n);

int mulMat(const double** mat1, const double** mat2, double** matRes, int rows1, int cols1rows2, int cols2);

int calcEquationSolution(const double** matA, const double** matB, double** matX, int n);

int calcMatTransformation(const double** src, const double** dst, double** matT, int n);

#endif


