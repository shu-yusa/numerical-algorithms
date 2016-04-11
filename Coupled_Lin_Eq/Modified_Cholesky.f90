!====================================================================!
!     Title  : Modified_Cholesky.f90                                 !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-4-25-Fri                                         !
!====================================================================!
      program  Modified_Cholesly
      implicit none
      integer :: i, j, N
      integer, allocatable :: mm(:)
      real(8), allocatable :: a(:,:), b(:), D(:)
      real(8), allocatable :: x(:), y(:)

      call  Title

      open(10, file='Matrix_4x4_SO_2', action='read')
      read(10,*) N
      allocate( a(N,N), b(N), x(N), y(N), D(N), mm(N) )
      do i=1, N
          read(10,*) (a(i,j), j=1, N), b(i)
      end do
      call  Matrix(N, a, b)

      do i=1, N
          D(i) = a(i, i)
      end do
      x = b

      call  LDLT_Decomp(N, a, b, mm, D, x)
      
      close(10)
      call  Solution(N, mm, x)

      stop
      end program
!====================================================================!
      subroutine LDLT_Decomp(N, L, b, mm, D, x)
      implicit none
      integer :: i, k
      integer, intent(in) :: N
      integer, intent(out) :: mm(N)
      real(8), intent(in) :: b(N)
      real(8), intent(inout) :: L(N,N), D(N)
      real(8), intent(out) :: x(N)
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

      do k=1, N
        x(k) = b(mm(k)) - dot_product(L(k,1:k-1), x(1:k-1))
      end do
      do k=N, 1, -1
        x(k) = x(k) / D(k) - dot_product(L(k+1:N,k), x(k+1:N))
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
      subroutine  Title
      
      write(6,*)
      write(6,*) '**************************'
      write(6,*) ' Modified Cholesky Method '
      write(6,*) '**************************'
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  Matrix(N, a, b)
      implicit none
      integer, intent(in) :: N
      integer :: i, j
      real(8), intent(in) :: a(N,N), b(N)
      character(len=30) :: FM='(x,a,4f6.2,a,i1,a,f7.2,a)'

      do i=1,N
        write(6,FM)  '|', (a(i, j), j=1, N), &
                   & ' || x(',i,') | = |',b(i), ' |'
      end do
      
      return
      end subroutine
!====================================================================!
      subroutine  Solution(N, mm, x)
      implicit none
      integer, intent(in) :: N, mm(N)
      integer :: i, t(N)
      real(8), intent(in) :: x(N)
      character(len=20) :: FM='(x,a,i1,a,f7.2)'

      do i=1, N
        t(mm(i)) = i
      end do

      write(6,*)
      do i=1, N
        write(6,FM) 'x(',i,') =',x(t(i))
      end do
      write(6,*)

      return
      end subroutine
