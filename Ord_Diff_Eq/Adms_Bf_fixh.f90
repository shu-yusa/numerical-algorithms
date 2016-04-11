!====================================================================!
!     Title  : Adms_Bf_fixh.f90                                      !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-3-24-Mon                                         !
!     Last modified : 2008-3-26-Wed                                  !
!                                                                    !
!     A program which solves ordinary differential equations by      !
!     using Adams-Bashforth's formula.                               !
!     The step width is fixed and is calculated from Fehlberg's      !
!     formula. In this formula, we can adjust the step number 'N'    !
!     depending on the value of local discretization error.          !
!                                                                    !
!     ** Structure of this program **                                !
!                                                                    !
!     module      Com_var                                            !
!     program     Adams_bashforth                                    !
!     subroutine  show_1           contained -> Adams_bashforth      !
!     subroutine  show_2           contained -> Adams_bashforth      !
!     subroutine  Init_Steps                                         !
!     subroutine  Cal_gamma                                          !
!     subroutine  Back_Diff                                          !
!     subroutine  Make_x_y_PHIe                                      !
!     subroutine  Next_B_Diff                                        !
!     subroutine  Estim_Error                                        !
!     subroutine  Manip_N                                            !
!     subroutine  funct                                              !
!     module      Fehlb_Constants                                    !
!     subroutine  Fberg                                              !
!     subroutine  func                                               !
!     function    Solu                                               !
!     subroutine  line                                               !
!     subroutine  Show_Eq                                            !
!     subroutine  Title                                              !
!====================================================================!
      module Com_var
!--------------------------------------------------------------------!
!     Definition global constants.                                   !
!--------------------------------------------------------------------!
      implicit none
      integer, parameter :: NMAX  = 15, NVAR = 1
      integer, parameter :: ALPHA = 101
      real*8  :: epsa, epsr
      parameter(epsa=1.0d-300, epsr=1.0d-15)
      end module
      
!====================================================================!
      program Adams_bashforth
!--------------------------------------------------------------------!
!    Main program.                                                   !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer :: k, N = 5, step = 0
      real*8  :: x, y(NVAR), h, f(NVAR), Err(-1:2), xout = 0.0d0
      real*8  :: PHI(NVAR, 0:NMAX+1), gam(0:NMAX+1), Error
      real*8, external :: Solu
      character ::  FM1*50
      parameter(FM1='(1x,f9.6,3x,i5,3x, 1pd14.7,3x,1pd13.5,3x,i2)')

      call Title
      call Show_Eq

      x = 0.0d0
      y = 1.0d0
      h = 1.0d0

      call  Init_Steps(N, x, y, h, PHI)
      call  show_1
      call  Cal_gamma(gam)
      call  Back_Diff(N, x, y, PHI)
      Err(1) = 0.0d0

      do while(x < 1.0d0) 
          call  Make_x_y_PHIe(N, x, y, h, gam, PHI)
          call  funct(1, x, y ,f)
          call  Next_B_Diff(N, f, PHI)
          call  Estim_Error(N, PHI, gam, Err, h)
          call  show_2
          call  Manip_N(N, Err, y)
      end do

      contains
!====================================================================!
      subroutine  show_1
!--------------------------------------------------------------------!
!     A subroutine which prints out the results calculated from      !
!     Fehlberg's formula.                                            !
!--------------------------------------------------------------------!
      character :: FM2*50
      parameter(FM2='(5x,"x",8x,"step",10x,"y",13x,"Error", 8x,"N")')

      write(6,FM2)
      call line
      x = - h
      do k=0, N
          x = x + h
          if (x >= xout) then
              Error = PHI(1,k) - Solu(1, x)
              write(6,FM1) x, step, PHI(1,k), Error, N
              xout = xout + 0.1d0
          end if
      end do

      end  subroutine
!====================================================================!
      subroutine  show_2
!--------------------------------------------------------------------!
!     A subroutine which prints out the results calculated from      !
!     Adams-bashforth's formula.                                     !
!--------------------------------------------------------------------!
          step = step + 1
      if (x >= xout) then
          Error = y(1) - Solu(1, x)
          write(6,FM1) x, step, y(1), Error, N
          xout = xout + 0.1d0
          step = 0
      end if

      end  subroutine
      end  program
     
!====================================================================!
      subroutine  Init_Steps(N, x, y, h, PHI)
!--------------------------------------------------------------------!
!     We calculate the first 'N' steps by Fehlberg's formula.        !
!     To calculate 'y_n', we have to know 'y_0' ~ 'y_n-1'. Therefore !
!     we use Fehlberg's formula.                                     !
!     After calling this subroutine, we have :                       !
!       x = x_N = x0 + N * h                                         !
!       y = y_N                                                      !
!       PHI(:, 0) = y_0                                              !
!       PHI(:, 1) = y_1                                              !
!       ...                                                          !
!       PHI(:, N) = y_N                                              !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer, intent(in) :: N
      real*8, intent(inout) :: x, y(NVAR), h
      real*8, intent(out)   :: PHI(NVAR, 0:NMAX+1)
      integer :: ini_step = 0
      real*8 :: xmin, h0

      xmin = x
      h0 = h
      PHI(:, 0) = y(:)

 A:   do while (ini_step < N)
          if (ini_step == 0) then
              x    = xmin
              y(:) = PHI(:, 0)
          end if
          PHI(:, ini_step) = y(:)

          call  Fberg(x, y, h, 0.1d0)
          
          h = min(h, 0.1d0)
          ini_step = ini_step + 1
          if (h < h0) then            ! If the new step width is smaller 
              ini_step = 0            ! than previous one, recalculate 
              h0 = h                  ! from first step.
          else 
              h = h0
          end if
      end do   A

      PHI(:, N) = y(:)

      write(6,'(7x,a,1pd14.7)') 'step size = ', h
      write(6,*)
      
      return
      end  subroutine

!====================================================================!
      subroutine Cal_gamma(gam)
!--------------------------------------------------------------------!
!     Calculation of gamma_k (0 <= k <= NMAX+1)                      !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer :: k, L
      real*8, intent(out) :: gam(0:NMAX+1)

      gam(0) = 1.0d0
      do k=1, NMAX+1
          gam(k) = 1.0d0
          do L=1, k
              gam(k) = gam(k) - gam(k - L) / (L + 1)
          end do
      end do
      
      return
      end subroutine

!====================================================================!
      subroutine Back_Diff(N, x, y, PHI)
!--------------------------------------------------------------------!
!     A subroutine which makes Back Difference table.                !
!     We compute difference between PHIs and swap them               !
!     After the calculation, we have :                               !
!                                                                    !
!        PHI(:, 0) = f_N = f(x_N, y_N)                               !
!        PHI(:, 1) = nabla    f_N                                    !
!        PHI(:, 2) = nabla^2  f_N                                    !
!        ...                                                         !
!        PHI(: ,N) = nabla^N  f_N                                    !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer :: j, k, L, N
      real*8  :: x, xx, y, h, w(NVAR), f(NVAR)
      real*8, intent(inout) :: PHI(NVAR, 0:NMAX+1)

      xx = x - N * h   ! = xmin

      do k=0, N
          do j=1, NVAR
              call  funct(j, xx, PHI(:,k), f)
              PHI(j, k) = f(j)
          end do
          xx = xx + h
      end do

      do L=1, N
          do k=0, N-L
              PHI(:, k) = PHI(:, k+1) - PHI(:, k)
          end do
      end do

      L = 0
      k = N
      do while(L < k)
          w(:)  = PHI(:, L)
          PHI(:, L) = PHI(:, k)
          PHI(:, k) = w(:)
          L = L + 1
          k = k - 1
      end do
      
      return
      end  subroutine

!====================================================================!
      subroutine  Make_x_y_PHIe(N, x, y, h, gam, PHI)
!--------------------------------------------------------------------!
!     Calculation of 'x_n+1' and 'y_n+1' from 'x_n' and 'y_n'.       !
!     We also make                                                   !
!                                                                    !
!     PHI(:, k) = nabla^k f_N + nabla^(k+1) f_N + ... + nabla^N f_N  !
!     (k = 0 ~ N-1)                                                  !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer :: k, N
      real*8, intent(inout) :: x, y(NVAR), PHI(NVAR, 0:NMAX+1)
      real*8, intent(in)    :: h, gam(0:NMAX+1)
      real*8  :: x0, S(NVAR)

      x0 = x
      x  = x + h
      S = 0.0d0

      S(:) = gam(N-1) * PHI(:, N-1)
      do k=N-2, 0, -1
          S(:)  = S(:) + gam(k) * PHI(:, k)
          PHI(:, k) = PHI(:, k+1) + PHI(:, k)
      end do
      y(:) = y(:) + h * S(:)

      return
      end  subroutine

!====================================================================!
      subroutine  Next_B_Diff(N, f, PHI)
!--------------------------------------------------------------------!
!     This subroutine calculates the next back difference which is   !
!     needed for calculating the next step.                          !
!     After the calculation, we have :                               !
!                                                                    !
!      PHI(:, 0)   = f_(n+1)                                         !
!      ...                                                           !
!      PHI(:, k)   = nabla^k     f_(n+1)                             !
!      ...                                                           !
!      PHI(:, N)   = nabla^N     f_(n+1)                             !
!      PHI(:, N+1) = nabla^(N+1) f_(n+1)                             !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer, intent(in)   :: N
      real*8, intent(inout) :: f(NVAR), PHI(NVAR, 0:NMAX+1)
      integer :: k

      PHI(:, N+1) = f(:) - PHI(:, 0) - PHI(:, N)
      PHI(:, N)   = f(:) - PHI(:, 0)
      PHI(:, 0)   = f(:)

      do k=1, N-1
          PHI(:, k) = PHI(:, k) + PHI(:, N)
      end do

      return
      end  subroutine

!====================================================================!
      subroutine  Estim_Error(N, PHI, gam, Err, h)
!--------------------------------------------------------------------!
!     Estimation of Error :                                          !
!      e_(N+2), e_(N+1), e_N, e_(N-1) .                              !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer, intent(in) :: N
      real*8, intent(in)  :: h, gam(0:NMAX+1), PHI(NVAR, 0:NMAX+1)
      real*8, intent(out) :: Err(-1:2)
      
      if (N > 2) Err(-1) = maxval( abs(PHI(:, N-2)) ) * h * gam(N-2)
      if (N > 1) Err(0)  = maxval( abs(PHI(:, N-1)) ) * h * gam(N-1)
      Err(1) = maxval( abs(PHI(:, N))   ) * h * gam(N)
      Err(2) = maxval( abs(PHI(:, N+1)) ) * h * gam(N+1)

      return
      end  subroutine 

!====================================================================!
      subroutine  Manip_N(N, Err, y)
!--------------------------------------------------------------------!
!     Manipulation of N according to the value of errors.            !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer, intent(inout) :: N
      real*8, intent(in)     :: y(NVAR)
      real*8, intent(inout)  :: Err(-1:2)
      integer :: Nnew
      real*8  :: eps

      eps = epsa + epsr * maxval( abs(y) )

      if (Err(1) > eps)  then 
          write(6,*) 'crash!'
          write(6,*) 'Err =', Err(1)
          write(6,*) 'eps =', eps
          stop
      else 
          Nnew = N
          if (N > 1) then
              if (N == 2) then
                  if ( Err(0) <= 0.5*Err(1) ) then
                      Nnew = N - 1
                      Err(1) = Err(0)
                  end if
              else
                  if ( max(Err(0), Err(-1)) < Err(1) ) then
                      Nnew = N - 1
                      Err(1) = Err(0)
                  end if
              end if
         end if

         if ( (N < NMAX).and.(N == Nnew) ) then
             if ( N == 1 ) then
                 if ( Err(2) <= 0.5*Err(1) ) then
                     Nnew = N + 1
                     Err(1) = Err(0)
                 end if
             else
                 if ( Err(2) <= min(Err(1), Err(0)) ) then
                     Nnew = N
                     Err(1) = Err(0)
                 end if
             end if
         end if

         N = Nnew
      end if
      
      return
      end  subroutine

!====================================================================!
      subroutine  funct(i, x, y, f)
!--------------------------------------------------------------------!
!     Definition of function.                                        !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer :: i
      real*8  :: x, y(NVAR), f(NVAR)

      if (i > NVAR) then
          write(6, *) 'Error in the number of function'
          write(6,*) 'i    =',i
          write(6,*) 'NVAR =',NVAR
          stop
      end if

      select case(i)
          case(1);  f(1) = - ALPHA * y(1)
      end select

      return
      end subroutine

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
      real*8, public  :: w(1:6), dw(1:6), alph(2:6), beta(2:6, 5)
      real*8, private ::  w1,  w2,  w3,  w4,  w5,  w6
      real*8, private :: dw1, dw2, dw3, dw4, dw5, dw6
      real*8, private :: alph2,  alph3,  alph4,  alph5,  alph6
      real*8, private :: beta21, beta22, beta23, beta24, beta25,     &
                       & beta31, beta32, beta33, beta34, beta35,     &
                       & beta41, beta42, beta43, beta44, beta45,     &
                       & beta51, beta52, beta53, beta54, beta55,     &
                       & beta61, beta62, beta63, beta64, beta65
      parameter(                                                     &
        & w1 =    16.d0/135.0d0, w2 =     0.0d0,                     &
        & w3 = 6656.0d0/12825.0d0, w4 = 28561.0d0/56430.0d0,         &         
        & w5 = - 0.18d0,         w6 =     2.0d0/   55.0d0 )
      parameter(                                                     &
        & dw1 =     1.0d0/ 360.0d0, dw2 =    0.0d0,                  &          
        & dw3 = - 128.0d0/4275.0d0, dw4 = - 2197.0d0/75240.0d0,      &
        & dw5 =    0.02d0,          dw6 =    2.0d0/   55.0d0 )
      parameter(                                                     &
        & alph2 = 0.25,          alph3 = 0.375d0,                    &
        & alph4 = 12.0d0/13.0d0, alph5 = 1.0d0,                      &
        & alph6 = 0.5d0 )
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
      subroutine  Fberg(x, y, h, xout)
!--------------------------------------------------------------------!
!     A subroutine of Fehlberg's formula.                            !
!     We compute next step width 'h'.                                !
!     Since the step width is not fixed, we change 'h' to 'xout - x' !
!     if 'h' plus 'x' is greater than 'xout'. Then we can print out  !
!     the results at 'xout'.                                         !
!--------------------------------------------------------------------!
      use  Fehlb_Constants
      use  Com_var
      implicit none
      integer :: i, j, k
      real*8, intent(in)    :: xout
      real*8, intent(inout) :: x, y(NVAR), h
      real*8 :: x0, y0(NVAR), kk(NVAR, 1:6)
      real*8 :: yerr, dy(NVAR), ynorm, eps
      real*8, external :: f

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
          ynorm = max(ynorm, abs(y(i)))
          yerr  = max(yerr,  abs(dy(i)))
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
      use Com_var
      implicit none
      integer :: i
      real*8  :: x, y(*), f(*)

      f(1) = - ALPHA * y(1)
      !f(2) = - 10 * ALPHA * y(2)

      return
      end subroutine
      
!====================================================================!
      real*8  function Solu(i, x)
!--------------------------------------------------------------------!
!     Exact solution of the differential equation.                   !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer :: i
      real*8  :: x
      
      select case(i)
          case(1);  Solu = dexp(- ALPHA * x)
          case(2);  Solu = dexp(- 10 * ALPHA * x)
      end select

      return 
      end
 
!====================================================================!
      subroutine  line
      
      write(6,*) '-----------------------------------------------&
                 &----------'
      
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
      write(6,*) '         *****************************'
      write(6,*) '           Adams-Bashforth (h,fixed)  '
      write(6,*) '         *****************************'

      return 
      end subroutine
