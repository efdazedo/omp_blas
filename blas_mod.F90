      module blas_mod
      implicit none
      interface
          subroutine daxpy(n,a,x,incx,y,incy)
          implicit none
          integer :: n, incx,incy
          real*8 :: a, x(*), y(*)
          end subroutine daxpy

          real*8 function ddot(n,x,incx,y,incy)
          implicit none
          integer :: n,incx,incy
          real*8 :: x(*),y(*)
          end function ddot

          real*8 function dnrm2(n,x,incx)
          implicit none
          integer :: n,incx
          real*8 :: x(*)
          end function dnrm2

        subroutine dgemv(trans,m,n,alpha,A,lda,x,incx,beta,y,incy)
        implicit none
        character :: trans
        integer :: m,n,lda,incx,incy
        real*8 :: alpha, beta, A(lda,*),x(*),y(*)
        end subroutine dgemv

        subroutine dgemm(transA,transB,m,n,k,                            &
     &   alpha,A,lda,B,ldb,beta,C,ldc )
        implicit none
        character transA,transB
        integer m,n,k,lda,ldb,ldc
        real*8 alpha,beta
        real*8 A(lda,*),B(ldb,*),C(ldc,*)
        end subroutine dgemm

        end interface

        end module blas_mod

