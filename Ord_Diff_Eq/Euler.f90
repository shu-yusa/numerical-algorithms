!====================================================================!
!     Title  : Euler.f90                                             !
!     Author : Yusa Shusaku                                          !
!     Data   : 2008-3-10-Mon                                         !
!====================================================================!
      program Euler
      implicit none
      integer :: i, n
      real*8,allocatable  :: x(:), y(:)
      real*8  :: a, b, h
      real*8,external :: f
      character(len=30) :: F100
      parameter(F100='(f19.15,f19.15)')
      parameter(a=0.0d0, b=1.0d0, h=0.1d0)

      n = NINT( (b - a) / h )
      allocate( x(1:n+1), y(1:n+1))
      x = 0.0d0
      y = 1.0d0

      write(6,F100) x(1) , y(1)
      do i=1, n
         x(i+1) = x(i) + h
         y(i+1) = y(i) + h*f(x(i), y(i))
         write(6,F100) x(i+1), y(i+1)
      end do
      !call pict
      
      end program

!====================================================================!
      real*8  function f(x,y)
      implicit none
      real*8,intent(in) :: x, y
      
      f = - x * x * y 

      return 
      end function

!====================================================================!
      subroutine pict

      write(6,*)
      write(6,*) '|************************************|'
      write(6,*) '|                                    |'
      write(6,*) '|           dy                       |'
      write(6,*) '|          ---- = - x^2 * y          |'
      write(6,*) '|           dx                       |'
      write(6,*) '|                                    |'
      write(6,*) '|          y(0) = 1                  |'
      write(6,*) '|                                    |'
      write(6,*) '**************************************'
      write(6,*)

      end subroutine
