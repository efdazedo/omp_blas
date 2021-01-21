      subroutine dgemm_omp(transA,transB, mm,nn,kk,                      &
     &  alpha, A,lda, B,ldb, beta, C, ldc )
!$omp declare target
      implicit none
      character, intent(in) :: transA, transB
      integer, intent(in) :: mm,nn,kk
      integer, intent(in) :: lda,ldb,ldc
      real*8 :: alpha, beta
      real*8, intent(in) :: A(lda,*), B(ldb,*)
      real*8, intent(inout) :: C(ldc,*)

      real*8, parameter :: zero = 0.0d0
      integer, parameter :: nb = 64

      real*8 :: cij, aik, bkj
      integer :: istart,iend, jstart,jend
      integer :: i,j,k
      logical :: is_transA, is_transB
      logical :: is_beta_zero
      integer :: m,n

      m = mm
      n = nn
      is_transA = (transA.eq.'t').or.(transA.eq.'T')
      is_transB = (transB.eq.'t').or.(transB.eq.'T')

      is_beta_zero = (beta.eq.zero)

!$omp parallel do simd collapse(2) private(i,j)
      do j=1,n
      do i=1,m
        if (is_beta_zero) then
                C(i,j) = zero
        else
                C(i,j) = beta * C(i,j)
        endif
      enddo
      enddo

      do jstart=1,n,nb
      do istart=1,m,nb
       jend = min(n,jstart+nb-1)
       iend = min(m,istart+nb-1)

!$omp  parallel do simd collapse(2)                                       &
!$omp& private(i,j,k,   cij,aik,bkj)
      do j=jstart,jend
      do i=istart,iend
        cij = zero
        do k=1,kk
          aik = merge( A(k,i), A(i,k), is_transA )
          bkj = merge( B(j,k), B(k,j), is_transB )
          cij = cij + aik * bkj
        enddo

        C(i,j) = C(i,j) + alpha * cij

       enddo
       enddo

       enddo
       enddo

       return
       end subroutine dgemm_omp
