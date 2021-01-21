      subroutine dcopy_omp( n, x, incx, y, incy )
!$omp declare target
      implicit none
      integer, intent(in) :: n, incx, incy
      real*8, intent(in) :: x(*)
      real*8, intent(inout) :: y(*)

      integer :: i, ix, iy

      if ((incx.eq.1).and.(incy.eq.1)) then
!$omp parallel do simd 
        do i=1,n
         y(i) = x(i)
        enddo
      else
!$omp parallel do simd private(ix,iy)
        do i=1,n
          ix = 1 + (i-1)*incx
          iy = 1 + (i-1)*incy
          y(iy) = x(ix)
        enddo
      endif

      return
      end subroutine dcopy_omp
             
