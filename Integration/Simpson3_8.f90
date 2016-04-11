!====================================================================!
!     Title  : Simpson3_8.f90                                        !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-5-25-Sun                                         !
!====================================================================!
      module  Global_Constant
      real*8, parameter :: epsa=1.0d-200, epsr=1.0d-15
      end module
!====================================================================!
      program main
      implicit none
      real*8 :: a, b, S, dS
      real*8 :: PI=3.141592653589793d0
            
      a = - 1.0d0
      b = 1.0d0
      call Simpson3_8(a, b, S, dS)
      write(6,*)
      write(6,*) 'S   =', S
      write(6,*) 'Err =', abs(S - PI)
      write(6,*)

      end program
!====================================================================!
      subroutine  Simpson3_8(a, b, S, dS)
      use Global_Constant, only : epsa, epsr
      implicit none
      integer :: i, j, jmax, N, MXHLF
      real*8 :: a, b, h
      real*8 :: S, S0, Sw, dS
      real*8, external :: f
      parameter(MXHLF=100)

      h = b - a
      S = h / 16.0d0 *(f(a+0.5d0*h) + f(a+h) &
     &    + 3.0d0 * (f(a+2.0d0*h/3.0d0) + f(a+5.0d0*h/6.0d0)))
      N = 1.0d0
      jmax = 1
      i = 0
      h = 0.5d0 * h

      do 
         i = i + 1
         jmax = 2 * jmax + 1
         h = 0.5d0 * h
         
         Sw = 0.0d0
         do j=1, jmax-2, 2 
            Sw = Sw + (f(a+dble(j)*h) + f(a+dble(j+2)*h) &
        &        + 3.0d0 * (f(a+dble(3*j+2)*h/3.0d0) &
        &                 + f(a+dble(3*j+4)*h/3.0d0)))
         end do
         Sw = Sw + (f(a+dble(jmax)*h) &
        &     + 3.0d0 * f(a+dble(3*jmax+2)*h/3.0d0))

         S0 = S
         S = 0.5d0 * S + 0.125d0 * h * Sw
         dS = abs(S - S0)
         if (dS < epsa + epsr * abs(S)) then
            return
         end if
         if(i > MXHLF) then
              write(6,*) 'The caculation did not converge.'
              write(6,*) 'S =', S
              exit
         end if
         if (h < 5.0d-7) exit
      end do


      return
      end subroutine
!====================================================================!
      function f(x) result(testf)
      implicit none
      real*8,intent(in) :: x
      real*8 :: testf

      testf = 2.0d0 / (1.0d0 + x * x)

      return
      end function
