!====================================================================!
!     Title  : LU_decomp.f90                                         !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-4-24-Thu                                         !
!     Last modified : 2008-4-24-Thu                                  !
!====================================================================!
      program main
      implicit none
      integer :: i, j, N
      real*8, allocatable :: a(:,:), b(:)
      real*8, allocatable :: x(:)

      call  Title

      open(10, file='Matrix_4x4_SO_2', action='read')
      read(10,*) N
      allocate( a(N,N), b(N), x(N) )
      do i=1, N
          read(10,*) (a(i,j), j=1, N), b(i)
      end do

      call  Matrix(N, a, b)
      call  LU_decomp(N, a, b, x)
      call  Results(N, a, b, x)

      close(10)

      stop
      end program
!====================================================================!
      subroutine  Matrix(N, a, b)
      implicit none
      integer, intent(in) :: N
      integer :: i
      real*8, intent(in) :: a(N,N), b(N)
      character(len=30) :: FM1='(x,a,4f7.2,a,i1,a,f7.2,a)'

      do i=1, N
          write(6,FM1) '|', a(i,:), ' || x(',i,') | = |', b(i), ' |'
      end do
      write(6,*)
      
      return
      end subroutine
!====================================================================!
      subroutine  Results(N, a, b, x)
      implicit none
      integer, intent(in) :: N
      integer :: i
      real*8, intent(in) :: a(N,N), b(N), x(N)
      character(len=30) :: FM='(x, a, i1, a, f7.2)'

      do i=1, N
          write(6,FM) 'x(', i, ') =',x(i)
      end do
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  Title

      write(6,*)
      write(6,*) '*******************'
      write(6,*) '  LU decompositon  '
      write(6,*) '*******************'
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  LU_decomp(N, a, b, x)
      implicit none
      integer, intent(in) :: N 
      integer :: mm(N)
      real*8, intent(in) :: a(N,N), b(N)
      real*8, intent(inout) :: x(N)
      real*8 :: p
            
      call  Decompose(N, a, p, mm)
      call  For_subst(N, a, b, mm, x)
      call  Back_subst(N, a, x)

      return
      end subroutine
!====================================================================!
      subroutine  Decompose(N, a, p, mm)
      implicit none
      integer :: i, j, k
      integer, intent(in) :: N
      integer, intent(out) :: mm(N)
      real*8, intent(inout) :: a(N,N)
      real*8, intent(in) :: p

      do i=1, N
          mm(i) = i
      end do
      
      do k=1, N-1
          call Par_Pivot(N, a, k, mm, p)
          a(k, k+1:N) = a(k, k+1:N) / p
          do i=k+1, N
              a(i, k+1:N) = a(i, k+1:N) - a(i, k) * a(k, k+1:N)
          end do
      end do

      return
      end subroutine
!====================================================================!
      subroutine  Par_Pivot(N, a, k, mm, p)
      implicit none
      integer, intent(in) :: k, N
      integer :: iw, j, L
      real*8, intent(out) :: p
      integer, intent(inout) :: mm(N)
      real*8, intent(inout) :: a(N, N)
      real*8 :: w(N)

      L = k
      p = a(k, k)

      do j=k+1, N
          if ( abs(p) < abs(a(j, k)) ) then
              L = j
              p = a(j, k)
          end if
      end do

      if (L /= k) then
          w(:)    = a(k, :)
          a(k, :) = a(L, :)
          a(L, :) = w(:)

          iw    = mm(k)
          mm(k) = mm(L)
          mm(L) = iw
      end if

      return
      end subroutine
!====================================================================!
      subroutine For_subst(N, a, b, mm, x)
      implicit none
      integer :: j, k
      integer, intent(in) :: mm(N), N
      real*8, intent(out) :: x(N)
      real*8, intent(in) :: a(N, N), b(N)
      
      do k=1, N
          x(k) = ( b(mm(k)) - DOT_PRODUCT(a(k,1:k-1), x(1:k-1)) ) &
                  / a(k,k)
      end do

      return
      end subroutine
!====================================================================!
      subroutine Back_subst(N, a, x)
      implicit none
      integer :: j, k
      integer, intent(in) :: N
      real*8, intent(in) :: a(N, N)
      real*8, intent(inout) :: x(N)

      do k=N, 1, -1
          x(k) = x(k) - Dot_Product(a(k,k+1:N), x(k+1:N))
      end do

      return
      end subroutine
