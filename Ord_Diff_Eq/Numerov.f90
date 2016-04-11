!====================================================================!
!     Title  : Numerov.f90                                           !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-6-19-Thu                                         !
!                                                                    !
!                                                                    !
!     *** Structure of this program ***                              !
!                                                                    !
!     module      Com_var                                            !
!     program     main                                               !
!     subroutine  Next_y                                             !
!     subroutine  Fix_h                                              !
!     subroutine  Lin_Correct                                        !
!     function    S                                                  !
!     function    f                                                  !
!     subroutine  Show_Result                                        !
!     subroutine  Title                                              !
!     subroutine  Show_Eq                                            !
!     module      Cash_Karp_Constants                                !
!     subroutine  Cash_Karp                                          !
!     subroutine  func                                               !
!                                                                    !
!====================================================================!
      module Com_var
      implicit none
      integer, parameter :: NVAR=2
      real(8) :: epsa, epsr
      real(8), parameter :: PI=3.1415926535897932d0
      parameter(epsa=1.0d-200, epsr=1.0d-15)
      end module
!====================================================================!
      program main
      implicit none
      integer :: i, n
      real(8) :: x, y(11), h, b, xout=2.0d0
      real(8), external :: PHI 

      call  Title
      call  Show_Eq
      x = 0.0d0
      !call  Fix_h(x, h)
      h = 0.05d0
      b = 20.0d0
      y(3) = 0.0d0 
      y(2) = PHI(h) 
      x = x + 2.0d0 * h

      write(6,*)
      do 
         call  Next_y(x, y(1:3), h)
         call  Show_Result(x, y(1), h, xout)
         x = x + h
         if (x > b) exit
         y = CSHIFT(y,-1)
      end do
      write(6,*)

      stop
      end program
!====================================================================!
      subroutine  Next_y(x, y, h)
      implicit none
      integer :: i, j, k, n
      real(8)  :: wh, rhs, ynm2, ynm1, denom
      real(8), intent(in) :: x, h
      real(8), intent(inout) :: y(3)
      real(8), external :: f, S

      wh = h * h / 12.0d0
      rhs = wh * (S(x) + 10.0d0 * S(x-h) + S(x-2.0d0*h))
      ynm1 = (2.0d0 - 10.0d0 * wh * f(x-h) ** 2) * y(2)
      ynm2 = (1.0d0 + wh * f(x-2.0d0*h) ** 2) * y(3)
      y(1) = (rhs + ynm1 - ynm2) * (1.0d0 - wh * f(x) ** 2)

      return
      end subroutine
!====================================================================!
      subroutine  Fix_h(x, h)
      implicit none
      real(8), intent(in) :: x
      real(8), intent(inout) :: h
      real(8) :: x0, y(2)

      x0 = x
      y(1) = 0.0d0
      y(2) = 0.5d0
      h = 0.1d0
      call  Cash_Karp(x0, y, h, 0.1d0)

      return
      end subroutine
!====================================================================!
      subroutine  Lin_Correct(x, y, h)
      implicit none
      integer :: i
      real(8), intent(in) :: x, h
      real(8), intent(inout) :: y(11)
      real(8) :: a, b

      a = 0.1d0 * (y(1) - y(11)) / h
      do i=1, 11
         y(i) = y(i) - a * (x + (1.0d0 - dble(i)) * h)
      end do

      return
      end subroutine
!====================================================================!
      function S(x)  result(T)
      implicit none
      real(8), intent(in) :: x
      real(8) :: T

      T = - 0.5d0 * x * exp(- x)

      end function
!====================================================================!
      function f(x)  result(g)
      implicit none
      real(8), intent(in) :: x
      real(8) :: g
      
      g = 0.0d0

      end function
!====================================================================!
      function PHI(x)  result(S)
      implicit none
      real(8), intent(in) :: x
      real(8) :: S
      
      S = 1.0d0 - 0.5d0 * (x + 2.0d0) * exp(-x)

      end function
!====================================================================!
      subroutine  Show_Result(x, y, h, xout)
      implicit none
      real(8), intent(in) :: x, y, h
      real(8), intent(inout) :: xout
      real(8), external :: PHI
      real(8) :: err
      character(len=30) :: FM='(x,a,f10.7,a,f9.6,a,1pd15.8)'
      
      if (x + 0.0001d0 >= xout) then
         err = abs(y - PHI(x))
         write(6,FM) 'x =',x, ' | y =', y, &
     &               ' | err =',err
         xout = xout + 2.0d0
      end if

      return
      end subroutine
!====================================================================!
      subroutine  Title

      write(6,*) 
      write(6,*) '             ****************'
      write(6,*) '              Numerov Method '
      write(6,*) '             ****************'
      write(6,*) 

      return
      end subroutine
!====================================================================!
      subroutine  Show_Eq

      write(6,*) '  |******************************************|'
      write(6,*) '  | Differential Equation :                  |'
      write(6,*) '  |                                          |'
      write(6,*) '  |  d^2y     　 1                           |'
      write(6,*) '  | ------  = - --- * x * exp(-x)            |'
      write(6,*) '  |  dx^2  　    2                           |'
      write(6,*) '  |                                          |'
      write(6,*) '  | Conditions :                             |'
      write(6,*) '  |              dy                          |'
      write(6,*) '  |   y(0) = 0, ----(0) = 0.5                |'
      write(6,*) '  |              dx                          |'
      write(6,*) '  |                                          |'
      write(6,*) '  |  Exact Solution :                        |'
      write(6,*) '  |               1                          |'
      write(6,*) '  |   y(x) = 1 - --- * (x + 2) * exp(-x)     |'
      write(6,*) '  |               2                          |'
      write(6,*) '  |                                          |'
      write(6,*) '  |******************************************|'
      write(6,*)

      return
      end subroutine
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
      real(8), public  :: w(1:6), dw(1:6), alph(2:6), beta(2:6, 5)
      real(8), private ::  w1,  w2,  w3,  w4,  w5,  w6
      real(8), private :: dw1, dw2, dw3, dw4, dw5, dw6
      real(8), private :: alph2,  alph3,  alph4,  alph5,  alph6
      real(8), private :: beta21, beta22, beta23, beta24, beta25,     &
                       & beta31, beta32, beta33, beta34, beta35,     &
                       & beta41, beta42, beta43, beta44, beta45,     &
                       & beta51, beta52, beta53, beta54, beta55,     &
                       & beta61, beta62, beta63, beta64, beta65
      parameter(                                                     &
        & w1 =  37.0d0/378.0d0, w2 =   0.0d0,                        &
        & w3 = 250.0d0/621.0d0, w4 = 125.0d0/ 594.0d0,               &
        & w5 =   0.0d0,         w6 = 512.0d0/1771.0d0 )
      parameter(                                                     &
        & dw1 = - 277.0d0/ 64512.0d0, dw2 =      0.0d0,              &
        & dw3 =  6925.0d0/370944.0d0, dw4 = - 6925.0d0/202752.0d0,   &
        & dw5 = - 277.0d0/ 14336.0d0, dw6 =    277.0d0/  7084.0d0 )
      parameter(                                                     &
        & alph2 = 0.2d0,   alph3 = 0.3d0,                            &
        & alph4 = 0.6d0,   alph5 = 1.0d0,                            &
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
      subroutine  Cash_Karp(x, y, h, xout)
!--------------------------------------------------------------------!
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
      real(8) :: x0, y0(NVAR), kk(NVAR, 1:6)
      real(8) :: yerr, dy(NVAR), ynorm, eps

      if ( x + h > xout )  h = xout - x
      x0 = x
      y0 = y
      
      call func(x, y, kk(:,1))

  AA: do j=2, 6
          x = x0 + alph(j) * h
   B:     do i=1, NVAR
              y(i) = 0.0d0
              do k=1, j-1
                  y(i) = y(i) + beta(j,k) * kk(i,k)
              end do
              y(i) = y0(i) + h * y(i)
          end do   B

          call func(x, y, kk(:,j))
      
      end do   AA

      x = x0 + h
      ynorm = 0.0d0
      yerr  = 0.0d0
  CC: do i=1, NVAR
          y(i)  = 0.0d0
          dy(i) = 0.0d0
   D:     do k=1, 6
              y(i)  = y(i)  +  w(k) * kk(i,k)
              dy(i) = dy(i) + dw(k) * kk(i,k)
          end do   D
          y(i)  = y0(i) + h * y(i)
          ynorm = max(ynorm, abs(y(i)))
          yerr  = max(yerr,  abs(dy(i)))
      end do   CC

      eps  = epsa + epsr * ynorm
      h    = 0.9d0 * h * (eps / (h * yerr) ) ** 0.2d0
      
      return
      end  subroutine
!===================================================================!
      subroutine func(x, y, g)
!--------------------------------------------------------------------!
!     Definition of the RHS of differential equations.               !
!     T = 0 corresponds to the potential without LS coupling.        !
!     T = 1 corresponds to the potential with LS coupling coupled    !
!     in parallel.                                                   !
!     T = 2 corresponds to the potential with LS coupling coupled    !
!     in antiparallel.                                               !
!--------------------------------------------------------------------!
      use Com_var, only : NVAR
      implicit none
      integer :: i, L
      real(8), intent(inout) :: x, y(NVAR), g(NVAR)
      real(8) :: U
      real(8), external :: S, f

      g(1) = y(2)
      g(2) = S(x) - f(x) ** 2 * y(1)

      return
      end subroutine
