      module blas_omp
      implicit none
      contains
#include "daxpy_omp.F90"
#include "ddot_omp.F90"
#include "dgemv_omp.F90"
      end module blas_omp
