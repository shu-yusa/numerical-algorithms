!====================================================================!
!     Title  : Taylor.f90                                            !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-3-12-Wed                                         !
!                                                                    !
!     A program which solves an ordinary differential equation by    !
!     using 'Taylor expansion method'.                               !
!====================================================================!
      program Taylor
      implicit none
      integer :: i, n, m, MMAX
      real*8  :: c(0:3), h, a, b
      real*8, external :: f
      real*8, allocatable :: x(:), y(:)
      parameter(a=0.0d0, b=1.0d0, h=0.1d0)
      parameter(MMAX=4)

      n = NINT((b - a)/h)
      allocate( x(1:n+1), y(1:n+1) )
      x(1) = 0.0d0
      y(1) = 1.0d0
      
      do i=1,n
          x(i+1) = x(i) + h
          c(0) = y(i)
          c(1) = - x(i)*x(i)*c(0)
          c(2) = - 0.5d0 * x(i) * ( 2*c(0) + x(i)*c(1) )
          y(i+1) = c(0) + h*( c(1) + h*c(2) ) 

!  Coefficients for m>=3 .
          do m=3,MMAX
              c(3) = - ( c(0) + x(i)*( 2*c(1) + x(i)*c(2) ) ) / m
              y(i+1) = y(i+1) + c(3)*h**m
              c(0) = c(1) 
              c(1) = c(2)
              c(2) = c(3)
          end do
      
      end do

      call Display(x,y,n)

      stop
      end program
             
!====================================================================!
      real*8  function f(x,y)
      implicit none
      real*8 :: x, y

      f = - x * x * y

      return
      end function

!====================================================================!
      subroutine  Display(x,y,n)
      implicit none
      integer, intent(in) :: n
      real*8,  intent(in) :: x(1:n+1), y(1:n+1)
      integer :: i
      character(len=20) :: FM1, FM2
      parameter(FM1='(f19.15,f19.15)')
      parameter(FM2='(1x,a,8x,a,18x,a)')

      write(6,*) 
      write(6,*) '#**********************************#'
      write(6,*) '# Differential Equ.                #'
      write(6,*) '#                                  #'
      write(6,*) '#         dy                       #'
      write(6,*) '#        ---- = - x^2*y            #'
      write(6,*) '#         dx                       #'
      write(6,*) '#                                  #'
      write(6,*) '#**********************************#'

      write(6,*)
      write(6,'(5x,a)') '#***********#'
      write(6,'(5x,a)') '#  Results  #'
      write(6,'(5x,a)') '#***********#'
      write(6,*)
 
      write(6,FM2) '#','x','y'
      do i=1,n+1
          write(6,FM1) x(i), y(i)  
      end do

      return
      end subroutine
