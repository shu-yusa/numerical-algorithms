!====================================================================!
!     Title  : Inv_mat_Cholesky.f90                                  !
!     Author : Yusa Shusaku                                          !
!     Date   : 2009-8-9-Sun                                          !
!====================================================================!
      program main
      implicit none
      integer :: i, j, N
      real(8), allocatable :: a(:,:), b(:), a_inv(:,:), c(:,:)
      real(8), allocatable :: x(:)
      character(30), parameter :: fm='(1x,a,i5,a,1pd22.15)'

      N = 500
      allocate( a(N,N), b(N), x(N), c(N,N), a_inv(N,N) )
      do i=1, N
        do j=1, i
          a(i,j) = dble(N - i + 1)
        end do
        do j=i+1, N
          a(i,j) = dble(N - j + 1)
        end do
      end do
      c = a

      if (N <= 9) then
        write(6,*) '** Original Matrix **'
        call Show_Matrix(N, a)
      end if
      call Inv_Mat(N, a, a_inv)
      if (N <= 9) then
        write(6,*) '** Inverse Matrix **'
        call Show_Matrix(N, a_inv)
        write(6,*) 'Check the product of them'
      end if
      a = matmul(c, a_inv)
      if (N <= 9) then
        call Show_Matrix(N, a)
      end if
      write(6,*) '** Check Diagonal element **'
      do i=1, N
        write(6,fm) 'i =',i,':', a(i,i)
      end do

      stop
      end program
!====================================================================!
      subroutine Inv_mat(N, a, a_inv)
      implicit none
      integer, intent(in) :: N
      integer :: i, k, mm(N)
      real(8), intent(in) :: a(N,N)
      real(8), intent(out) :: a_inv(N,N)
      real(8), dimension(N) :: b, x, D
      real(8), allocatable :: L(:,:)

      allocate(L(N,N))
      L = a
      forall (i=1:N) D(i) = a(i,i)
      call LDLT_Decomp(N, L, D, mm)
      b = 0.0d0
      do i=1, N
        b(i) = 1.0d0
        do k=1, N
          x(k) = b(mm(k)) - dot_product(L(k,1:k-1), x(1:k-1))
        end do
        do k=N, 1, -1
          x(k) = x(k) / D(k) - dot_product(L(k+1:N,k), x(k+1:N))
        end do
        a_inv(:,i) = x(:)
        b(i) = 0.0d0
      end do
      deallocate(L)

      return
      end subroutine
!====================================================================!
      subroutine LDLT_Decomp(N, L, D, mm)
      implicit none
      integer :: i, k
      integer, intent(in) :: N
      integer, intent(out) :: mm(N)
      real(8), intent(inout) :: L(N,N), D(N)
      real(8) :: u(N)
      
      forall (i=1:N) mm(i) = i

      do k=1, N
        call  Pivotting(N, k, L, mm, D(k))
        do i=1, k-1
          u(i) = L(k,i) - dot_product(u(1:i-1), L(i,1:i-1))
          L(k,i) = u(i) / D(i)
        end do
        D(k) = D(k) - dot_product(u(1:k-1), L(k,1:k-1))
      end do

      return
      end subroutine
!====================================================================!
      subroutine  Pivotting(N, k, a, mm, p)
      implicit none
      integer, intent(in) :: N
      integer :: j, L, iw
      integer, intent(in) :: k
      integer, intent(out) :: mm(N)
      real(8), intent(inout) :: a(N,N)
      real(8), intent(out) :: p
      real(8) :: w(N)

      L = k
      p = a(k,k)

      do j=k+1, N
        if (abs(p) < abs(a(j,j))) then
          L = j
          p = a(j,j)
        end if
      end do

      if (L /= k) then
        w(:)   = a(k,:)
        a(k,:) = a(L,:)
        a(L,:) = w(:)
            
        w(:)   = a(:,k)
        a(:,k) = a(:,L)
        a(:,L) = w(:)

        iw    = mm(k)
        mm(k) = mm(L)
        mm(L) = iw
      end if

      return
      end subroutine
!====================================================================!
      subroutine  Show_Matrix(N, a)
      implicit none
      integer, intent(in) :: N
      integer :: i, j
      real(8), intent(in) :: a(N,N)
      character(len=30) :: FM

      FM = '(x,a,' // CHAR(48+N) // 'f7.3,a)'
      do i=1, N
          write(6,FM) ' |', (a(i,j), j=1, N), ' |'
      end do
      write(6,*)

      return
      end subroutine
