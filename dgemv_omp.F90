      subroutine dgemv_omp(trans,m,n,alpha,A,lda,x,incx,beta,y,incy)
      implicit none
      character, intent(in) ::  trans
      integer, intent(in) :: m,n,lda,incx,incy
      real*8, intent(in) :: alpha, beta, A(lda,*), x(*)
      real*8, intent(inout) :: y(*)

      logical :: is_trans
      integer :: i,j,ix,iy, nrowx,nrowy
      real*8 :: aij, Ax

      is_trans = (trans.eq.'T').or.(trans.eq.'t').or.                    &
     &           (trans.eq.'C').or.(trans.eq.'c')

      
      nrowy = merge( n, m, is_trans )
      nrowx = merge( m, n, is_trans )

!$omp parallel do simd private(i,iy)
      do i=1,nrowy
        iy = 1 + (i-1)*incy
        if (beta.eq.0) then
                y(iy) = 0
        else
                y(iy) = beta*y(iy)
        endif
      enddo


!      -------------------------------------
!      y(:) += alpha * op(A(1:m,1:n)) * x(:)
!      -------------------------------------
!$omp    parallel do simd private(i,j,ix,iy,Ax,aij)
          do i=1,nrowy

           Ax = 0
           do j=1,nrowx
             aij  = merge( A(j,i), A(i,j), is_trans )
             ix = 1 + (j-1)*incx
             Ax = Ax + aij*x(ix)
           enddo

           iy = 1 + (i-1)*incy
           y(iy) = y(iy) + alpha*Ax
          enddo

          return
          end subroutine dgemv_omp
