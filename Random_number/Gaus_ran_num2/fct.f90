!======================================================================!
      function f(x) result(g)
!----------------------------------------------------------------------!
!     Definition of Gaussian function which will be integrated.        !
!----------------------------------------------------------------------!
      implicit none
      real*8, intent(in) :: x
      real*8, parameter :: PI=3.141592653589793d0
      real*8 :: g

      g = exp(- 0.5d0 * x * x) / sqrt(2.0d0 * PI)

      return
      end function
