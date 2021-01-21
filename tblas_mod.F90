       program tblas_mod
#ifdef _OPENMP
       use blas_omp
#else
       use blas_mod
#define ddot_omp ddot
#define daxpy_omp daxpy
#define dgemv_omp dgemv
#endif

       implicit none
       integer, parameter :: n = 100*1000
       integer, parameter :: ncase = 80*10
       real*8 :: x(n,ncase), y(n,ncase) 
       real*8 :: xy(ncase)
       real*8 :: alpha
       integer :: i,icase, incx,incy

       integer :: tstart,tend,count_rate
       real*8 :: elapsed_time

!$omp target data map(from:x,y,xy)

!$omp  target teams
!$omp  distribute parallel do simd collapse(2)
       do icase=1,ncase
       do i=1,n
         x(i,icase) = dble(i + icase)/dble(n*ncase)
         y(i,icase) = dble(-2*i + icase)/dble(n*ncase)
       enddo
       enddo
!$omp end target teams

       incx = 1
       incy = 1
       alpha = 1.1

!$omp taskwait
!$omp barrier
       call system_clock(tstart,count_rate)

!$omp target teams
!$omp distribute 
       do icase=1,ncase

        call daxpy_omp(n, alpha, x(:,icase), incx, y(:,icase), incy)
        xy(icase) = ddot_omp(n,x(:,icase),incx, y(:,icase),incy)
       enddo
!$omp end target teams

!$omp taskwait
!$omp barrier
       call system_clock(tend,count_rate)
       elapsed_time = dble(tend-tstart)/dble(count_rate)
       print*,'time is ', elapsed_time

!$omp end target data 
       print*,'sum(x) ', sum(x)
       print*,'sum(y) ', sum(y)
       print*,'sum(xy) ', sum(xy)

       stop
       end program tblas_mod


