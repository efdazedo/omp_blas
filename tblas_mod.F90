       program tblas_mod
#ifdef _OPENMP
       use blas_omp,                                                     &
     &   ddot => ddot_omp,                                               &
     &   dnrm2 => dnrm2_omp,                                             &
     &   dcopy => dcopy_omp,                                             &
!     &   dgemv => dgemv_omp,                                             &
!     &   dgemm => dgemm_omp,
     &   daxpy => daxpy_omp

#define dgemv dgemv_omp
#define dgemm dgemm_omp
#else
       use blas_mod
#endif

       implicit none
       integer, parameter :: n = 1024*2
       integer, parameter :: ncase = n
       real*8 :: x(n,ncase), y(n,ncase), z(n,ncase)
       real*8 :: xy(ncase)
       real*8 :: ynorm(ncase), znorm(ncase)
       real*8 :: alpha, beta
       integer :: i,icase, inc1,inc2
       integer :: istart,iend
       integer :: mm,nn,kk,ld1,ld2,ld3
       character :: transA, transB
       integer, parameter :: nb = 64

       integer :: tstart,tend,count_rate
       real*8 :: elapsed_time

!$omp target data map(from:x,y,z,xy,ynorm,znorm)

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
       beta = 1.2
       transA = 'N'

!$omp taskwait
!$omp barrier
       call system_clock(tstart,count_rate)

!$omp target teams
!$omp distribute  private(mm,nn,ld1)
       do icase=1,ncase
        call dcopy(n, z(:,icase), inc1, x(:,icase), inc2 )
        call daxpy(n, alpha, x(:,icase), inc1, y(:,icase), inc2)

        mm = size(z,1)
        nn = size(z,2)
        ld1 = size(z,1)
        call dgemv( transA, mm,nn,alpha, z,ld1,                          &
     &     x(:,icase),inc1,beta,y(:,icase),inc2 )

        xy(icase) = ddot(n,x(:,icase),inc1, y(:,icase),inc2)
        ynorm(icase) = dnrm2(n,y(:,icase),inc1)
       enddo
!$omp end target teams

        transA = 'N'
        transB = 'N'
        ld1 = size(x,1)
        ld2 = size(y,1)
        ld3 = size(z,1)
!$omp target teams
!$omp distribute private(istart,iend,mm,nn,kk)
      do istart=1,ncase,nb
        iend = min(ncase,istart+nb-1)
        mm = size(z,1)
        nn = iend-istart+1
        kk = size(x,2)
        call  dgemm(transA,transB,mm,nn,kk,                              &
     &          alpha, x,ld1, y,ld2,                                     &
     &          beta,  z(:,istart),ld3)
       enddo
!$omp end target teams

!$omp target teams
!$omp distribute 
       do icase=1,ncase
        znorm(icase) = dnrm2( size(z,1), z(:,icase), inc1)
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
       print*,'sum(ynorm) ',sum(ynorm)
       print*,'sum(znorm) ',sum(znorm)

       stop
       end program tblas_mod


