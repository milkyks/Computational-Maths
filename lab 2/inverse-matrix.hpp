#ifndef INVERSE_MATRIX_HPP
#define INVERSE_MATRIX_HPP

#include <math.h>

void solve(int n, double **a, double *b, int *ipvt);

void decomp(int n, double **a, double *cond, int *ipvt, double *work);

#endif
