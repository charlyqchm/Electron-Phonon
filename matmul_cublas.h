#ifndef MATMUL_CUBLAS
#define MATMUL_CUBLAS

#include "atom.h"
#include <cuComplex.h>
#include <cublas_v2.h>
#include <ctime>
#include <cstdlib>
void matcublas(complex<double> *matA, complex<double> *matB,
               complex<double> *matC, UNINT dim);

#endif
