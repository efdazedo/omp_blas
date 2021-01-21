       program tblas_mod
#ifdef _OPENMP
       use blas_omp,                                                     &
     &   ddot => ddot_omp,                                               &
     &   dcopy => dcopy_omp,                                             &
     &   dgemv => dgemv_omp,                                             &
     &   daxpy => daxpy_omp
#else
       use blas_mod
#endif

       implicit none
       integer, parameter :: n = 100*1000
       integer, parameter :: ncase = 80*10
       real*8 :: x(n,ncase), y(n,ncase), z(n,ncase)
       real*8 :: xy(ncase)
       real*8 :: alpha
       integer :: i,icase, inc1,inc2

       integer :: tstart,tend,count_rate
       real*8 :: elapsed_time

!$omp target data map(from:x,y,xy)

!$omp  target teams
!$omp  distribute parallel do simd collapse(2)
       do icase=1,ncase
       do i=1,n
         z(i,icase) = dble(i + icase)/dble(n*ncase)
         y(i,icase) = dble(-2*i + icase)/dble(n*ncase)
         x(i,icase) = 0
       enddo
       enddo
!$omp end target teams

       inc1 = 1
       inc2 = 1
       alpha = 1.1

!$omp taskwait
!$omp barrier
       call system_clock(tstart,count_rate)

!$omp target teams
!$omp distribute 
       do icase=1,ncase
        call dcopy(n, z(:,icase), inc1, x(:,icase), inc2 )
        call daxpy(n, alpha, x(:,icase), inc1, y(:,icase), inc2)
        xy(icase) = ddot(n,x(:,icase),inc1, y(:,icase),inc2)
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


