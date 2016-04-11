!====================================================================!
!     Title  : Jacobi.f90                                            !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-4-26-Sat                                         !
!     Last modified : 2008-4-27-Sun                                  !
!====================================================================!
      module  Global_Constant
      implicit none
      real*8, parameter :: epsr = 1.0d-15
      real*8, parameter :: epsa = 1.0d-100
      end module
!====================================================================!
      module  Guess_x
      contains
      subroutine  Guess_Solution(N, x)
      implicit none
      integer, intent(in) :: N
      real*8, intent(out) :: x(N)

      x(1) = -33.0d0
      x(2) = 9.0d0
      x(3) = 6.0d0

      return
      end subroutine
      end module
!====================================================================!
      program  Main
      implicit none
      integer :: i, j, N, Nmax, Step
      real*8, allocatable :: a(:,:), b(:), x(:)
      character(len=20) :: FNAME
      
      call Title

      FNAME = 'Matrix_4x4_SO'
      open(10, file=FNAME,action='read')
      read(10,*) N
      allocate(a(N,N), b(N), x(N))
      do i=1, N
          read(10,*) (a(i,j), j=1, N), b(i)
      end do

      call  Show_Matrix(N, a, b)
      call  Convergence(N, a, Nmax)
      call  Jacobi_2(N, Nmax, a, b, x, Step)
      call  Show_Result(N, a, b, x, Step)
      call  Show_Error(N, x)

      close(10)

      stop
      end program
!====================================================================!
      subroutine  Jacobi_1(N, Nmax, a, b, x, Step)
      use  Global_Constant
      use  Guess_x
      implicit none
      integer :: i, j, k
      integer, intent(in) :: N, Nmax
      integer, intent(out) :: Step
      real*8, intent(in) :: a(N,N)
      real*8, intent(out) :: x(N)
      real*8, intent(in) :: b(N)
      real*8 :: x0(N), dx, error, eps, r
      
      call  Guess_Solution(N, x)
      Step = 0
      do k=1, Nmax
          Step = Step + 1
          x0 = x
          error = 0.0d0
          do i=1, N
              r = b(i)
              do j=1, N
                  if (j == i) cycle
                  r = r - a(i,j) * x0(j)
              end do
              x(i) = r / a(i,i)
              dx = abs(x(i) - x0(i))
              eps = epsa + epsr * (abs(x0(i)) + abs(x(i)))
              error = max(error, dx / eps)
          end do
          if (error < 1.0d0) return
      end do

      stop 'Iteration did not converge.'

      return
      end subroutine
!====================================================================!
      subroutine  Jacobi_2(N, Nmax, a, b, x, Step)
      use  Global_Constant
      use  Guess_x
      implicit none
      integer, intent(in) :: N, Nmax
      integer :: i, j, k
      integer, intent(out) :: Step
      real*8, intent(inout) :: x(N)
      real*8, intent(in) :: a(N,N), b(N)
      real*8 :: error, dx, eps, r, x0(N)

      call  Guess_Solution(N, x)
      Step = 0
      do k=1, Nmax
          Step = Step + 1
          x0 = x
          error = 0.0d0
          do i=1, N
              r     = b(i) - Dot_Product(a(i,:), x0(:))
              dx    = r / a(i,i)
              x(i)  = x0(i) + dx
              eps   = epsa + epsr * (abs(x(i)) + abs(x0(i)))
              error = max(error, abs(dx) / eps)
          end do
          if (error < 1.0d0)  return
      end do
      
      stop 'Iteration did not converge.'

      return
      end subroutine
!====================================================================!
      subroutine  Convergence(N, a, Nmax)
      use Global_Constant
      implicit none
      integer :: i, j
      integer, intent(in) :: N
      integer, intent(out) :: Nmax
      real*8 :: anorm, amax
      real*8, intent(in) :: a(N,N)
      logical :: Converge
      character(len=20) :: FM='(x,a,i4)'
      character(len=20) :: FM2='(x,a,f7.4)'

      anorm = 0.0d0
      do i=1, N
          amax = 0.0d0
          do j=1, N
              if (j == i) cycle
              amax = amax + abs(a(i, j))
          end do
          anorm = max(anorm, amax / abs(a(i,i)))
      end do

      Converge = (anorm < 1.0d0)
      if ( Converge ) then
          write(6,*) 'Converge :', Converge
          write(6,FM2) 'anorm =',anorm
      else 
          write(6,*) 'Converge :', Converge
          write(6,FM2) 'anorm =',anorm
          write(6,*) 'Iteration may not converges.'
      end if

      Nmax = int(log(epsr) / log(anorm))
      write(6,FM) 'Nmax =', Nmax
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  Title

      write(6,*)
      write(6,*) '*************************'
      write(6,*) ' Jacobi Iteration Method '
      write(6,*) '*************************'
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  Show_Result(N, a, b, x, Step)
      implicit none
      integer, intent(in) :: N, Step
      integer :: i, j
      real*8, intent(in) :: a(N,N), b(N), x(N)
      character(len=30) :: FM2='(x,a,i1,a,f11.7)'
      character(len=30) :: FM ='(x,a,i3,a)'

      write(6,FM) 'Iteration :', Step, ' times'
      write(6,*)
      do i=1, N
          write(6,FM2) 'x(',i,') =', x(i)
      end do
      write(6,*) 

      return
      end subroutine
!====================================================================!
      subroutine  Show_Matrix(N, a, b)
      implicit none
      integer, intent(in) :: N
      integer :: i, j
      real*8, intent(in) :: a(N,N), b(n)
      character(len=30) :: FM='(x,a,4f7.2,a,i1,a,f7.2,a)'
      
      do i=1, N
          write(6,FM) '|', (a(i,j), j=1, N), &
                     &' || x(',i,') | = |', b(i), ' |'
      end do
      write(6,*) 

      return
      end subroutine
!====================================================================!
      subroutine  Show_Error(N, x)
      implicit none
      integer, intent(in) :: N
      real*8, intent(in) :: x(N)
      character(len=30) :: FM = '(x,a,1pd12.5)'

      write(6,FM) 'Err(1) = ', abs(x(1) - 1.0d0)
      write(6,FM) 'Err(2) = ', abs(x(2) - 2.0d0)
      write(6,FM) 'Err(3) = ', abs(x(3) + 1.0d0)
      write(6,FM) 'Err(4) = ', abs(x(4) - 1.0d0)
      write(6,*)

      return
      end subroutine

