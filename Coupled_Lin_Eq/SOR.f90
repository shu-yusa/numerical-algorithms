!====================================================================!
!     Title  : SOR.f90                                               !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-4-27-Sun                                         !
!     Last modified : 2008-4-27-Sun                                  !
!====================================================================!
      module  Global_Constant
      implicit none
      real*8, parameter :: epsr = 1.0d-15
      real*8, parameter :: epsa = 1.0d-100
      end module
!====================================================================!
      subroutine  Guess_Solution(N, Nmax, a, x, w)
      implicit none
      integer, intent(in) :: N
      integer, intent(out) :: Nmax
      integer :: i
      real*8, intent(in) :: a(N,N)
      real*8, intent(out) :: x(N), w
      real*8 :: anorm
      character(len=30) :: FM='(x,a,f8.5)',FM2='(x,a,i3)'

      do i=1, N
          x(i) = dble(i) * 2.0d0
      end do

      call  Convergence(N, a, anorm, Nmax)
      w = 2.0d0 / (1.0d0 + sqrt(1.0d0 - anorm * anorm))
      
      if (w > 0.0d0 .and. w < 2.0d0) then
          write(6,*) '|1 - omega| < 1'
      else
          write(6,*) '|1 - omega| > 1'
          write(6,*) 'Iteration will not converge.'
      end if

      write(6,FM) 'omega =', w
      write(6,FM2)  'Nmax  =', Nmax
      write(6,*)

      return
      end subroutine
!====================================================================!
      program  Main
      implicit none
      integer :: i, j, N, Nmax, Step
      real*8 :: anorm
      real*8, allocatable :: a(:,:), b(:), x(:)
      character(len=20) :: FNAME
      logical :: con

      FNAME = 'Matrix_4x4_SO'
      open(10, file=FNAME, action='read')
      read(10,*) N
      allocate(a(N,N), b(N), x(N))
      do i=1, N
          read(10,*) (a(i,j), j=1, N), b(i)
      end do

      call  Title
      call  Show_Matrix(N, a, b)
      call  SOR(N, Nmax, a, b, x, Step)
      call  Show_Result(N, a, b, x, Step)

      close(10)
      stop
      end program
!====================================================================!
      subroutine  SOR(N, Nmax, a, b, x, Step)
      use  Global_Constant
      implicit none
      integer, intent(in) :: N, Nmax
      integer, intent(out) :: Step
      integer :: i, j, k
      real*8, intent(inout) :: x(N)
      real*8, intent(in) :: a(N,N), b(N)
      real*8 :: error, dx, eps, r, xold, omega

      call  Guess_Solution(N, Nmax, a, x, omega)
      Step = 0
      do k=1, Nmax
          Step = Step + 1
          error = 0.0d0
          do i=1, N
              xold = x(i)
              r = b(i) - Dot_Product(a(i,1:N),x(1:N))
              dx = omega * r / a(i,i)
              x(i) = x(i) + dx
              eps = epsa + epsr * (abs(x(i)) + abs(xold))
              error = MAX(error, abs(dx) / eps)
          end do
          if (error < 1.0d0)  return
      end do
      
      stop 'Iteration did not converge.'

      return
      end subroutine
!====================================================================!
      subroutine  Convergence(N, a, anorm, Nmax)
      use Global_Constant
      implicit none
      integer :: i, j
      integer, intent(in) :: N
      integer, intent(out) :: Nmax
      real*8 :: amax
      real*8, intent(in) :: a(N,N)
      real*8, intent(out) :: anorm
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
      Nmax = int(log(epsr) / log(anorm))

      return
      end subroutine
!====================================================================!
      subroutine  Title

      write(6,*)
      write(6,*) '***********************************'
      write(6,*) ' Successive Over-Relaxation Method '
      write(6,*) '***********************************'
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  Show_Result(N, a, b, x, Step)
      implicit none
      integer, intent(in) :: N, Step
      integer :: i, j
      real*8, intent(in) :: a(N,N), b(N), x(N)
      character(len=30) :: FM2='(x,a,i1,a,f9.5)'
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
