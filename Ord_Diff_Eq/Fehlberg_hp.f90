!====================================================================!
!     Title  : Fehlberg.f90                                          !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-3-21-Fri                                         !
!     Last modified : 2008-3-22-Sat                                  !
!                                                                    !
!     A program which solves an ordinary differential equation by    !
!     using Fehlberg's formula.                                      !
!     In Fehlberg's formula, we can determine next step width such   !
!     that the error becomes smaller than the tollerable error.      !
!====================================================================!
      module Fehlb_Constants
!--------------------------------------------------------------------!
!     Definition of constants which appear in the Fehlberg's formula.!
!     Since substitution is not allowed in module subprogram and     !
!     operations such as division are not allowd in data sentence,   !
!     we have to define temporary variables to initialize arrays     !
!     used in other program units.                                   !
!--------------------------------------------------------------------!
      implicit none
      real(16), public  :: w(1:6), dw(1:6), alph(2:6), beta(2:6, 5)
      real(16), private ::  w1,  w2,  w3,  w4,  w5,  w6
      real(16), private :: dw1, dw2, dw3, dw4, dw5, dw6
      real(16), private :: alph2,  alph3,  alph4,  alph5,  alph6
      real(16), private :: beta21, beta22, beta23, beta24, beta25,     &
                       & beta31, beta32, beta33, beta34, beta35,     &
                       & beta41, beta42, beta43, beta44, beta45,     &
                       & beta51, beta52, beta53, beta54, beta55,     &
                       & beta61, beta62, beta63, beta64, beta65
      parameter(                                                     &
        & w1 =    16.d0/135.0d0, w2 =       0.0d0,                   &
        & w3 = 6656.0d0/12825.0d0, w4 = 28561.0d0/56430.0d0,         &         
        & w5 = - 0.18d0,         w6 =       2.0d0/   55.0d0 )
      parameter(                                                     &
        & dw1 =     1.0d0/ 360.0d0, dw2 =      0.0d0,                &
        & dw3 = - 128.0d0/4275.0d0, dw4 = - 2197.0d0/75240.0d0,      &
        & dw5 =    0.02d0,          dw6 =      2.0d0/   55.0d0 )
      parameter(                                                     &
        & alph2 = 0.25d0,        alph3 = 0.375d0,                    &
        & alph4 = 12.0d0/13.0d0, alph5 =   1.0d0,                    &
        & alph6 =  0.5d0 )
      parameter(                                                     &
        & beta21 =     0.25d0,          beta31 =     3.0d0/  32.0d0, &
        & beta32 =      9.0d0/  32.0d0, beta41 =  1932.0d0/2197.0d0, &      
        & beta42 = - 7200.0d0/2197.0d0, beta43 =  7296.0d0/2197.0d0, &
        & beta51 =    439.0d0/ 216.0d0, beta52 = -   8.0d0,          &
        & beta53 =   3680.0d0/ 513.0d0, beta54 = - 845.0d0/4104.0d0, &
        & beta61 = -    8.0d0/  27.0d0, beta62 =     2.0d0,          &
        & beta63 = - 3544.0d0/2565.0d0, beta64 =  1859.0d0/4104.0d0, &
        & beta65 = -  0.275d0 )
      data w    / w1,  w2, w3,  w4,  w5,  w6/
      data dw   /dw1, dw2,dw3, dw4, dw5, dw6/
      data alph /alph2, alph3, alph4, alph5, alph6/
      data beta(2,:) /beta21, 4*0.0d0/
      data beta(3,:) /beta31, beta32, 3*0.0d0/
      data beta(4,:) /beta41, beta42, beta43, 2*0.0d0/
      data beta(5,:) /beta51, beta52, beta53, beta54, 0.0d0/
      data beta(6,:) /beta61, beta62, beta63, beta64, beta65/

      end module Fehlb_Constants

!====================================================================!
      module  Com_var
!--------------------------------------------------------------------!
!     Definitions of parameters.                                     !
!--------------------------------------------------------------------!
      implicit none
      integer :: NVAR, ALPHA
      real(16)  :: epsa, epsr
      parameter(NVAR=1, epsa=1.0d-300, epsr=1.0d-15, ALPHA=10)
      end module Com_var

!====================================================================!
      program main
!--------------------------------------------------------------------!
!     Main program.                                                  !
!     We solve the equation in the range [0,1]. We print out the     !
!     result at the points 0.1*n with  n = 0, 1, ..., 10.            !
!     In Fehlberg's formula, we can adjust the step width so that    !
!     the error becomes small. However, an error in the initial step !
!     remains in subsequent steps. Therefore, we recalculate initial !
!     step after getting step width appropriate for it.              !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer :: step = 0
      real(16)  :: h, x, y(NVAR)
      real(16)  :: err, xout = 0.0d0
      real(16), external  :: Solu
      character(len=40) :: F1, F2
      parameter(F1='(3x,"x",4x,"STEP",9x,"y",14x,"ERROR")')
      parameter(F2='(1x,f4.1,3x,i4,3x,1pd15.8,3x,1pd13.5)')

      call Title
      call Show_Eq
      write(6,F1)  
      call line

!    Initial conditon
     
      x = 0.0d0
      y = 1.0d0
      h = 1.0d0

!    Determination of the initial step width.      
      
      call Fberg(x, y, h, 0.1d0)
      
      x = 0.0d0
      y = 1.0d0

 A:   do while(xout < 1.0d0)
          err = y(1) - Solu(1, x)
          
          if (x >= xout) then
              write(6, F2)  x, step, y(1), err
              xout = xout + 0.1d0
              step = 0
          end if
          
          call Fberg(x, y, h, xout)
         
          step = step + 1
      end do    A

      end program

!====================================================================!
      subroutine  Fberg(x, y, h, xout)
!--------------------------------------------------------------------!
!     A subroutine of Fehlberg's formula.                            !
!     We compute next step width 'h'.                                !
!     Since the step width is not fixed, we change 'h' to 'xout - x' !
!     if 'h' plus 'x' is greater than 'xout'. Then we can print out  !
!     the results at 'xout'.                                         !
!--------------------------------------------------------------------!
      use  Fehlb_Constants
      use  Com_var, only : epsa, epsr, NVAR
      implicit none
      integer :: i, j, k
      real(16), intent(in)    :: xout
      real(16), intent(inout) :: x, y(NVAR), h
      real(16) :: x0, y0(NVAR), kk(NVAR, 1:6)
      real(16) :: yerr, dy(NVAR), ynorm, eps
      !real(16), external :: f

      if ( x + h > xout )  h = xout - x
      x0 = x
      y0 = y
      
      call func(x, y, kk)

   A: do j=2, 6
          x = x0 + alph(j) * h
   B:     do i=1, NVAR
              y(i) = 0.0d0
              do k=1, j-1
                  y(i) = y(i) + beta(j,k) * kk(i,k)
              end do
              y(i) = y0(i) + h * y(i)
          end do   B

          call func(x, y, kk(:,j))
      
      end do   A

      x = x0 + h
      ynorm = 0.0d0
      yerr  = 0.0d0
   C: do i=1, NVAR
          y(i)  = 0.0d0
          dy(i) = 0.0d0
   D:     do k=1, 6
              y(i)  = y(i)  +  w(k) * kk(i,k)
              dy(i) = dy(i) + dw(k) * kk(i,k)
          end do   D
          y(i)  = y0(i) + h * y(i)
          ynorm = max(ynorm, Qabs(y(i)))
          yerr  = max(yerr,  Qabs(dy(i)))
      end do   C

      eps  = epsa + epsr * ynorm
      h    = 0.9d0 * h * (eps / (h * yerr) ) ** 0.2d0
      
      return
      end  subroutine

!====================================================================!
      subroutine func(x, y, f)
!--------------------------------------------------------------------!
!     Definition of function.                                        !
!--------------------------------------------------------------------!
      use Com_var, only : ALPHA, NVAR
      implicit none
      integer :: i
      real(16)  :: x, y(NVAR), f(NVAR)

      f(1) = - dble(ALPHA) * y(1)
    ! f(2) = - 10 * ALPHA * y(2)

      return
      end subroutine
      
!====================================================================!
      real(16)  function Solu(i, x)
!--------------------------------------------------------------------!
!     Exact solution of the differential equation.                   !
!--------------------------------------------------------------------!
      use Com_var, only : ALPHA
      implicit none
      integer :: i
      real(16)  :: x
      
      select case(i)
          case(1);  Solu = qexp(- dble(ALPHA) * x)
          case(2);  Solu = qexp(- 10.0d0 * dble(ALPHA) * x)
      end select

      return 
      end
 
!====================================================================!
      subroutine  line
      
      write(6,*) '-----------------------------------------------'
      
      return 
      end subroutine

!====================================================================!
      subroutine  Show_Eq
      use Com_var
      implicit none
      character :: Q*25
      parameter(Q='(1x,a,i3,8x,a)')

      write(6,*)
      write(6,*) '|**************************************|'
      write(6,*) '| Differential Equation :              |'
      write(6,*) '|                                      |'
      write(6,*) '|   dy                                 |'
      write(6,*) '|  ----  =  - a * x                    |'
      write(6,*) '|   dx                                 |'
      write(6,*) '|                                      |'   
      write(6,*) '| Conditions :                         |'
      write(6,*) '|                                      |'
      write(6,*) '|   0 <= x <= 1 ,  y(0) = 1            |'
      write(6,Q) '|   a = ', ALPHA ,'                    |'
      write(6,*) '|                                      |'
      write(6,*) '| Exact Solution :                     |'
      write(6,*) '|                                      |'
      write(6,*) '| y(x) = exp ( - a * x )               |'
      write(6,*) '|                                      |'
      write(6,*) '****************************************'
      write(6,*)

      end subroutine

!====================================================================!
      subroutine Title

      write(6,*)
      write(6,*) '         ********************'
      write(6,*) "          Fehlberg's formula "
      write(6,*) '         ********************'

      return 
      end subroutine
