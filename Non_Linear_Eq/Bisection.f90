      program main
      implicit none
      real*8 :: x
      x = 100.0d0

      call Bisection(x)
      write(6,*) 'x =', x

      stop
      end program
!======================================================================!
      subroutine Bisection(x)
      implicit none
      real*8, intent(inout) :: x
      real*8 :: xl, xr, x0
      real*8, parameter :: dx=2.0d0, epsr=1.0d-10
      real*8, external :: f
            
      do 
         if (f(x)*f(x+dx) > 0.0d0) then
            x = x + dx
         else
            xl = x  ;  xr = x + dx
            exit
      end if
      end do
         do 
            x0 = x
            x = 0.50d0 * (xl + xr)
            if (abs(x - x0) < epsr) exit
            if (f(x) > 0.0d0) then
               xr = x
            else 
               xl = x
            end if
         end do

      return         
      end subroutine
!======================================================================!
      function f(x)  result(g)
      implicit none
      real*8, intent(in) :: x
      real*8 :: g
      real*8, parameter :: mt=175.0d0, mb=5.0d0, mH=400.0d0

      g = sqrt(x*x + mt*mt) + sqrt(x*x + mb*mb) - mH

      return
      end function
