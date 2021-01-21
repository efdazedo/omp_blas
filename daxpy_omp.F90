      subroutine daxpy_omp( n, a, x,incx,y,incy)
!$omp declare target
      implicit none
      integer, intent(in) :: n, incx,incy
      real*8, intent(in) :: a
      real*8, intent(in) :: x(*)
      real*8, intent(inout) :: y(*)

      integer :: ix, iy, i

      if ((incx.eq.1).and.(incy.eq.1)) then
!$omp parallel do simd
              do i=1,n
                y(i) = a*x(i) + y(i)
              enddo

      else
!$omp parallel do simd private(ix,iy)
              do i=1,n
                ix = 1 + (i-1)*incx
                iy = 1 + (i-1)*incy
                y(iy) = a*x(ix) + y(iy)
              enddo
      endif
      return
      end subroutine daxpy_omp
