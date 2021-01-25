touch tblas_mod tblas_mod_nogpu
rm tblas_mod tblas_mod_nogpu
pgfortran  -o tblas_mod -fast -mp -gpu=cc70,cuda11.1 -Minfo=all  blas_omp.F90 tblas_mod.F90  -llapack -lblas
pgfortran  -o tblas_mod_nogpu  blas_mod.F90 tblas_mod.F90 -llapack -lblas
