      real*8 function dnrm2_omp(n,x,incx)
!$omp declare target
      implicit none
      integer, intent(in) :: n,incx
      real*8, intent(in) :: x(*)

      integer :: i,ix
      real*8 :: dsum

      dsum = 0
!$omp  parallel do simd                                                  &
!$omp& private(i,ix) reduction(+:dsum)
      do i=1,n
        ix = 1 + (i-1)*incx
        dsum = dsum + x(ix)*x(ix)
      enddo

      dnrm2_omp = sqrt( abs(dsum) )
      return
      end function dnrm2_omp
