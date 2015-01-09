#ifndef _EQUATION_SOLVE_H
#define _EQUATION_SOLVE_H


#define EPS 1e-100

int mallocMat(double*** matA, int rows, int cols);

int callocMat(double*** matA, int rows, int cols);

int calcMatDet(const double* const* mat, double* det, int n);

int revertMat(const double* const* mat, double** matRev, int rows, int cols);

int inverseMat(const double* const* mat, double** matInv, int n);

int mulMat(const double* const* mat1, const double* const* mat2, double** matRes, int rows1, int cols1rows2, int cols2);

int calcEquationSolution(const double* const* matA, const double* const* matB, double** matX, int n);

int calcMatTransformation(const double* const* src, const double* const* dst, double** matT, int n);

#endif


