#include "matmul_cublas.h"

void matcublas(complex<double> *matA, complex<double> *matB, complex<double> *matC,
            UNINT dim){

   int dim2 = dim * dim;
   // vector< complex<double> > aux_mat(dim2,0.0);
   cuDoubleComplex *dev_A, *dev_B, *dev_C;
   const cuDoubleComplex alf = make_cuDoubleComplex(1.0,0.0);
   const cuDoubleComplex bet = make_cuDoubleComplex(0.0, 0.0);
   const cuDoubleComplex *alpha = &alf;
   const cuDoubleComplex *beta = &bet;

   cudaMalloc((void**) &dev_A, dim2 * sizeof(cuDoubleComplex));
   cudaMalloc((void**) &dev_B, dim2 * sizeof(cuDoubleComplex));
   cudaMalloc((void**) &dev_C, dim2 * sizeof(cuDoubleComplex));

   cudaMemcpy(dev_A, matA, dim2 * sizeof(cuDoubleComplex),
              cudaMemcpyHostToDevice);
   cudaMemcpy(dev_B, matB, dim2 * sizeof(cuDoubleComplex),
              cudaMemcpyHostToDevice);
   cudaMemcpy(dev_C, matC, dim2 * sizeof(cuDoubleComplex),
              cudaMemcpyHostToDevice);
   // Create a handle for CUBLAS
  cublasHandle_t handle;
  cublasCreate(&handle);

  // Do the actual multiplication

  cublasZgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, dim, dim, dim, alpha, dev_A,
              dim, dev_B, dim, beta, dev_C, dim);

  // Destroy the handle
  cublasDestroy(handle);

  cudaMemcpy(matC, dev_C, dim2 * sizeof(cuDoubleComplex),
             cudaMemcpyDeviceToHost);

   //Free GPU memory
  cudaFree(dev_A);
  cudaFree(dev_B);
  cudaFree(dev_C);
}
