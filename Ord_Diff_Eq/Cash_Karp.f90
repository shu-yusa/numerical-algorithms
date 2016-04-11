!====================================================================!
!     Title  : Cash_Karp.f90                                         !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-3-22-Sat                                         !
!     Last modified : 2008-3-22-Sat                                  !
!                                                                    !
!     A program which solves an ordinary differential equation by    !
!     using Cash_Karp's formula.                                     !
!     In Cash_Karp's formula, we can determine next step width such  !
!     that the error becomes smaller than the tollerable error.      !
!====================================================================!
      module Cash_Karp_Constants
!--------------------------------------------------------------------!
!     Definition of constants which appear in the Cash_Karp's        !
!     formula.                                                       !
!     Since substitution is not allowed in module subprogram and     !
!     operations such as division are not allowd in data sentence,   !
!     we have to define temporary variables to initialize arrays     !
!     used in other program units.                                   !
!--------------------------------------------------------------------!
      implicit none
      real(8), public  :: w(1:6), dw(1:6), alph(2:6), beta(2:6,5)
      real(8), private ::  w1,  w2,  w3,  w4,  w5,  w6
      real(8), private :: dw1, dw2, dw3, dw4, dw5, dw6
      real(8), private :: alph2,  alph3,  alph4,  alph5,  alph6
      real(8), private :: beta21, beta22, beta23, beta24, beta25,    &
     &                    beta31, beta32, beta33, beta34, beta35,    &
     &                    beta41, beta42, beta43, beta44, beta45,    &
     &                    beta51, beta52, beta53, beta54, beta55,    &
     &                    beta61, beta62, beta63, beta64, beta65
      parameter(                                                     &
        & w1 =  37.0d0/378.0d0, w2 =   0.0d0,                        &
        & w3 = 250.0d0/621.0d0, w4 = 125.0d0/ 594.0d0,               &         
        & w5 =   0.0d0,         w6 = 512.0d0/1771.0d0 )
      parameter(                                                     &
        & dw1 = - 277.0d0/ 64512.0d0, dw2 =      0.0d0,              &
        & dw3 =  6925.0d0/370944.0d0, dw4 = - 6925.0d0/202752.0d0,   &
        & dw5 = - 277.0d0/ 14336.0d0, dw6 =    277.0d0/  7084.0d0 )
      parameter(                                                     &
        & alph2 =   0.2d0,  alph3 = 0.3d0,                           &
        & alph4 =   0.6d0,  alph5 = 1.0d0,                           &
        & alph6 = 0.875d0 )
      parameter(                                                     &
        & beta21 =    0.2d0,           beta31 = 0.075d0,             &
        & beta32 =  0.225d0,           beta41 =   0.3d0,             &      
        & beta42 = -  0.9d0,           beta43 =   1.2d0,             &
        & beta51 = - 11.0d0/   54.0d0, beta52 =   2.5d0,             &
        & beta53 = - 70.0d0/   27.0d0, beta54 =  35.0d0/27.0d0,      &
        & beta61 = 1631.0d0/55296.0d0, beta62 = 0.341796875d0,       &
        & beta63 =  575.0d0/13824.0d0, beta64 = 44275.0d0/110592.0d0,&
        & beta65 =  0.061767578125d0 )
      data w    / w1,  w2, w3,  w4,  w5,  w6/
      data dw   /dw1, dw2,dw3, dw4, dw5, dw6/
      data alph /alph2, alph3, alph4, alph5, alph6/
      data beta(2,:) /beta21, 4*0.0d0/
      data beta(3,:) /beta31, beta32, 3*0.0d0/
      data beta(4,:) /beta41, beta42, beta43, 2*0.0d0/
      data beta(5,:) /beta51, beta52, beta53, beta54, 0.0d0/
      data beta(6,:) /beta61, beta62, beta63, beta64, beta65/

      end module
!====================================================================!
      module  Com_var
!--------------------------------------------------------------------!
!     Definitions of parameters.                                     !
!--------------------------------------------------------------------!
      implicit none
      integer :: NVAR, ALPHA
      real(8)  :: epsa, epsr
      parameter(NVAR=1,epsa=1.0d-300,epsr=1.0d-15,ALPHA=10)
      end module
!====================================================================!
      program main
!--------------------------------------------------------------------!
!     Main program.                                                  !
!     We solve the equation in the range [0,1]. We print out the     !
!     result at the points 0.1*n with  n = 0, 1, ..., 10.            !
!     In Cash_Karp's formula, we can adjust the step width so that   !
!     the error becomes small. However, an error in the initial step !
!     remains in subsequent steps. Therefore, we recalculate initial !
!     step after getting step width appropriate for it.              !
!--------------------------------------------------------------------!
      use Cash_Karp_Constants
      use Com_var, only : NVAR
      implicit none
      integer :: step = 0
      real(8)  :: h, x, y(NVAR)
      real(8)  :: err, xout = 0.0d0
      real(8), external  :: Solu
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

!    Determination of initial step width.      
      
      call Next_y(x, y, h, 0.1d0)
      
      x = 0.0d0
      y = 1.0d0

 A:   do while(xout < 1.0d0)
          err = y(1) - Solu(x)
          if (x >= xout) then
              write(6, F2)  x, step, y(1), err
              xout = xout + 0.1d0
              step = 0
          end if
          call Next_y(x, y, h, xout)
          step = step + 1
      end do    A

      end program
!====================================================================!
      subroutine  Next_y(x, y, h, xout)
!--------------------------------------------------------------------!
!     A subroutine of Cash_Karp's formula.                           !
!     We compute next step width 'h'.                                !
!     Since the step width is not fixed, we change 'h' to 'xout - x' !
!     if 'h' plus 'x' is greater than 'xout'. Then we can print out  !
!     the results at 'xout'.                                         !
!--------------------------------------------------------------------!
      use  Cash_Karp_Constants
      use  Com_var, only : NVAR, epsa, epsr
      implicit none
      integer :: i, j, k
      real(8), intent(in)    :: xout
      real(8), intent(inout) :: x, y(NVAR), h
      real(8) :: x0, y0(NVAR), kk(NVAR, 1:6), yerr, dy, ynorm, eps

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

          call func(x, y, kk(1,j))
      
      end do   A

      x = x0 + h
      ynorm = 0.0d0
      yerr  = 0.0d0
   C: do i=1, NVAR
          y(i) = 0.0d0
          dy   = 0.0d0
   D:     do k=1, 6
              y(i) = y(i) +  w(k) * kk(i,k)
              dy   =  dy  + dw(k) * kk(i,k)
          end do   D
          y(i)  = y0(i) + h * y(i)
          ynorm = max(ynorm, abs(y(i)))
          yerr  = max(yerr,  abs(dy)  )
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
      real(8)  :: x, y(NVAR), f(NVAR)

      f(1) = - dble(ALPHA) * y(1)

      return
      end subroutine
!====================================================================!
      real(8)  function Solu(x)
!--------------------------------------------------------------------!
!     Exact solution of the differential equation.                   !
!--------------------------------------------------------------------!
      use Com_var, only : ALPHA
      implicit none
      real(8) :: x

      Solu = exp(- dble(ALPHA) * x)

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

      return
      end subroutine

!====================================================================!
      subroutine Title

      write(6,*)
      write(6,*) '         *********************'
      write(6,*) "          Cash-Karp's formula "
      write(6,*) '         *********************'

      return 
      end subroutine
