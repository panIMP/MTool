#pragma once

int mallocMat(double*** matA, int rows, int cols);

int calcMatDet(double** mat, double* det, int n);

int revertMat(double** mat, double** matRev, int rows, int cols);

int inverseMat(double** mat, double** matInv, int n);

int mulMat(double** mat1, double** mat2, double** matRes, int rows1, int cols1rows2, int cols2);

int calcEquationSolution(double** matA, double** matB, double** matX, int n);

int calcMatTransformation(double** src, double** dst, double** matT, int n);


