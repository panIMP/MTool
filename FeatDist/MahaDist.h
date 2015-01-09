#ifndef _MAHA_DIST_H_
#define _MAHA_DIST_H_


int calcMatCovar(const double* const* x, double** matCovar, int dim, int num);

int calcMahaDistance(const double* mat1, const double* mat2, const double* const* matCovarInv, int dim, double* dist);

int calcMahaDistance2(const double* mat1, const double* mat2, const double* const* matCovarInv, int dim, double* dist);

#endif