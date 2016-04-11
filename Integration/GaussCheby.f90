!====================================================================!
!     Title  : GaussCheby.f90                                        !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-5-25-Sun                                         !
!====================================================================!
      module Global_Constant
      implicit none
      integer, parameter :: MAXn=100
      real*8,  parameter :: epsa=1.0d-300, epsr=1.0d-15
      real*8,  parameter :: PI=3.141592653589793d0
      end module
!====================================================================!
      program Main
      use Global_Constant, only : epsr, epsa, MAXn, PI
      implicit none
      integer :: n, m
      real*8 :: S, S0
      real*8, external :: f, gam, Exact_Sol
      character(len=30) :: FM1='(x,a,1pd23.15)'

      m = 200
      S = PI * f(m,0.0d0)
      do n=2, MAXn
         S0 = S
         call  Integral(n, m, S)
         if (abs(S - S0) < epsa + epsr * abs(S)) exit
         if (n == MAXn) then
            write(6,FM1) 'S =', S
            stop 'Integration did not converge.'
         end if
      end do

      call  Title
      call  Show_Integral(m)
      call  Show_Result(n, m, S)

      end program
!====================================================================!
      subroutine Integral(n, m, S)
!--------------------------------------------------------------------!
!     Integration by Gauss-Chebyshev formula.                        !
!--------------------------------------------------------------------!
      use Global_Constant, only : PI
      implicit none
      integer, intent(in) :: n, m
      real*8, intent(out) :: S
      integer :: j
      real*8 :: a0
      real*8, external :: f

      a0 = 0.5d0 * PI / dble(n)

      S = 0.0d0
      do j=1, n
         S = S + f(m,cos(a0 * dble(2 * j - 1)))
      end do
      
      S = S * PI / dble(n)

      return
      end subroutine
!====================================================================!
      function f(m, x)  result(testf)
!--------------------------------------------------------------------!
!     Definition of test function whose integrated value we know.    !
!--------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: m
      real*8, intent(in)  :: x
      real*8 :: testf
      
      testf = x ** m
     
      return
      end function
!====================================================================!
      function  Exact_Sol(m)  result(S)
!--------------------------------------------------------------------!
!     Exact integrated value of test function.                       !
!--------------------------------------------------------------------!
      use Global_Constant, only : PI
      implicit none
      integer, intent(in) :: m
      real*8 :: S
      real*8, external :: gam

      S = sqrt(PI) * gam(dble(m+1)*0.5d0) / gam(dble(m+2)*0.5d0)

      return
      end function
!====================================================================!
      function gam(x)  result(g)
!--------------------------------------------------------------------!
!     Gamma function.                                                !
!     Argument can be an half-integer.                               !
!--------------------------------------------------------------------!
      use Global_Constant, only : PI
      implicit none
      real*8, intent(in) :: x
      real*8 :: g
      integer :: n, i
      logical :: c

      n = int(x)
      c = (x /= dble(n))
      if (c) then
         g = sqrt(PI)
         do i=1, n
            g = g * (x - dble(i))
         end do
      else 
         g = 1.0d0
         do i=2, n
            g = g * dble(i-1)
         end do
      end if

      return
      end function
!====================================================================!
      subroutine  Title

      write(6,*)
      write(6,*) '            *************************'
      write(6,*) '             Gauss-Chebyshev Formula '
      write(6,*) '            *************************'
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  Show_Integral(m)
      implicit none
      integer, intent(in) :: m
      real*8 ,external :: Exact_Sol
      character(len=30) :: P='(x,a,i3,a,37x,a)'
      character(len=30) :: Q='(x,a,1pd23.15,20x,a)'

      write(6,*) '|*************************************************|'
      write(6,*) '| Integral:                                       |'
      write(6,*) '|                                                 |'
      write(6,*) '|      1                                          |'
      write(6,*) '|     /                                           |'
      write(6,*) '|     [      x^m           sqrt(PI) * G((m+1/2))  |'
      write(6,*) '| S = | --------------- = ----------------------- |'
      write(6,*) '|     ]  sqrt(1 - x^2)        2 * G((m+2)/2)      |'
      write(6,*) '|     /                                           |' 
      write(6,*) '|     -1                                          |' 
      write(6,*) '|                                                 |'
      write(6,*) '| where G(x) is a gamma function.                 |'
      write(6,P) '| For m =',m,',',                                '|'
      write(6,*) '|                                                 |'
      write(6,Q) '|   S =', Exact_Sol(m),                          '|'
      write(6,*) '|                                                 |'
      write(6,*) '***************************************************'
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  Show_Result(n, m, S)
      implicit none
      integer, intent(in) :: n, m
      real*8, intent(in) :: S
      real*8, external :: Exact_Sol
      character(len=30) :: FM1='(2x,a,i3)'
      character(len=30) :: FM2='(3x,a,i3,a)', FM3='(3x,a,1pd23.15)'

      write(6,FM1) 'Computation with m = ',m
      write(6,*)
      write(6,*)   '*** Result ***'
      write(6,FM2) 'Iteration :', n ,' times'
      write(6,FM3) 'S   =', S
      write(6,FM3) 'ERR =', abs(S - Exact_Sol(m))
      write(6,*) 

      return
      end subroutine
