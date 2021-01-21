      real*8 function ddot_omp(n, x,incx,y,incy)
!$omp declare target
      implicit none
      integer, intent(in) :: n, incx, incy
      real*8, intent(in) :: x(*)
      real*8, intent(in) :: y(*)

      real*8 :: dsum
      integer :: i, ix,iy

      dsum = 0
!$omp parallel do simd private(ix,iy) reduction(+:dsum)
      do i=1,n
        ix = 1 + (i-1)*incx
        iy = 1 + (i-1)*incy
        dsum = dsum + y(iy) * x(ix)
      enddo

      ddot_omp = dsum
      return
      end function ddot_omp
