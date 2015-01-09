#ifndef _MAHA_DIST_H_
#define _MAHA_DIST_H_


int calcMatCovar(const double*** x, double** matCovar, int dim, int num);

int getMahaDistance(const double** mat1, const double** mat2, const double** matCovarInv, int dim, double* dist);

int getMahaDistance2(const double** mat1, const double** mat2, const double** matCovarInv, int dim, double* dist);

#endif