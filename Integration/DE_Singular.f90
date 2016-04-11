!====================================================================!
!     Title  : Dble_Exp.f90                                          !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-5-25-Sun                                         !
!====================================================================!
      module  Global_Constant
      implicit none
      integer, parameter :: MXHLF=100
      real*8, parameter :: epsa=1.0d-300, epsr=1.0d-15
      real*8, parameter :: PI=3.141592653589793d0
      real*8, parameter :: MAXH=6.0d0
      end module
!====================================================================!
      program Main
      use Global_Constant
      implicit none
      integer :: n
      real*8 :: a, b, S, dS

      a = 0.0d0
      b = 1.0d0
      call Title
      call DE(a, b, n, S, dS)
      call Show_Integral
      call Show_Result(n, S)

      end program
!====================================================================!
      subroutine  DE(a, b, n, S, dS)
      use Global_Constant, only : MAXH, MXHLF, epsa, epsr
      implicit none
      integer, intent(out) :: n
      integer :: nhlf, j, jmax
      real*8, intent(in) :: a, b
      real*8, intent(out) :: S, dS
      real*8 :: c1, c2, h, sh, exsh, chsh
      real*8 :: xj, xmj, Sw, fwj, S0, wj
      real*8, external :: f

      h = MAXH
      c1 = 0.5d0 * (b - a)
      c2 = 0.5d0 * (b + a)
      nhlf = 0
      jmax = 0
      sh   = sinh(MAXH)
      exsh = exp(-sh)
      chsh = cosh(sh)
      xj  = 0.5d0 * (a * exsh + b / exsh) / chsh
      xmj = 0.5d0 * (a / exsh + b * exsh) / chsh
      wj  = cosh(MAXH) * exsh * exsh / ((chsh * exsh) ** 2)
      S   = MAXH * (f(c2) + wj * (f(xj) + f(xmj)))

      do n=1, MXHLF 
         h = 0.5d0 * h 
         jmax = 2 * jmax + 1
         Sw = 0.0d0

         do j=jmax, 1, -2
            sh = sinh(dble(j)*h)
            exsh = exp(-sh)
            chsh = cosh(sh)
            xj  = 0.5d0 * (a * exsh + b / exsh) / chsh
            xmj = 0.5d0 * (a / exsh + b * exsh) / chsh
            fwj = cosh(dble(j)*h) * exsh * exsh / ((chsh * exsh)**2) &
     &             * (f(xj) + f(xmj))
            Sw = Sw + fwj
         end do   
         S0 = S
         S  = 0.5d0 * S + h * Sw
         dS = abs(S - S0)
         if (dS < epsa + epsr * abs(S)) then
            S = c1 * S
            return
         end if
      end do      

      write(6,*) 'S   =', S
      write(6,*) 'Err =', dS
      stop 'Integration did not converge.'

      end subroutine
!====================================================================!
      function  f(x)  result(testf)
      implicit none
      real*8, intent(in) :: x
      real*8 :: testf

      testf = 2.0d0 / sqrt(x * (2.0d0 - x))

      return
      end function
!====================================================================!
      subroutine  Title

      write(6,*)
      write(6,*) '         ****************************'
      write(6,*) '          Double-Exponential Formula '
      write(6,*) '         ****************************'
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  Show_Integral

      write(6,*) '|**********************************************|'
      write(6,*) '| Integral:                                    |'
      write(6,*) '|                                              |'
      write(6,*) '|      1                                       |'
      write(6,*) '|     /                                        |'
      write(6,*) '|     [         2                              |'
      write(6,*) '| S = | ----------------- = PI                 |'
      write(6,*) '|     ]   sqrt(x(1 - x))                       |'
      write(6,*) '|     /                                        |' 
      write(6,*) '|      0              = 3.14159265358797932... |'
      write(6,*) '|                                              |'
      write(6,*) '|                                              |'
      write(6,*) '************************************************'
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  Show_Result(n, S)
      use Global_Constant, only : PI
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: S
      character(len=30) :: FM2='(3x,a,i3,a)', FM3='(3x,a,1pd23.15)'

      write(6,*)
      write(6,*)   '*** Result ***'
      write(6,FM2) 'Iteration :', n ,' times'
      write(6,FM3) 'S   =', S
      write(6,FM3) 'ERR =', abs(S - PI)
      write(6,*) 

      return
      end subroutine
