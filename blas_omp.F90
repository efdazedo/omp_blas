      module blas_omp
      implicit none
      contains
#include "dcopy_omp.F90"
#include "dnrm2_omp.F90"
#include "daxpy_omp.F90"
#include "ddot_omp.F90"
#include "dgemv_omp.F90"
#include "dgemm_omp.F90"
      end module blas_omp
