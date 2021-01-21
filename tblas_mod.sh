touch tblas_mod tblas_mod_nogpu
rm tblas_mod tblas_mod_nogpu
ftn -o tblas_mod -homp blas_omp.F90 tblas_mod.F90 
ftn -o tblas_mod_nogpu  blas_mod.F90 tblas_mod.F90 
