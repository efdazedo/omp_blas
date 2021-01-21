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

        subroutine dgemv(trans,m,n,A,lda,alpha,x,incx,beta,y,incy)
        implicit none
        character :: trans
        integer :: m,n,lda,incx,incy
        real*8 :: alpha, beta, A(lda,*),x(*),y(*)
        end subroutine dgemv

        end interface

        end module blas_mod

