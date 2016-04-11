!======================================================================!
!     Title  : main.f90, fct.f90, GaussLegen_module.f90                !
!     Author : Yusa Shusaku                                            !
!     Date   : 2009-2-5-Thu                                            !
!     Last modified : 2009-2-5-Thu                                     !
!                                                                      !
!     This program generates the Gaussian-distributed random numbers   !
!     from random numbers distributed uniformally in the range [0,1].  !
!     The method is adopted from "Theory of nuclear reactions"         !
!     (by P.Froebrich & R.Lipperheide) p.398.                          !
!                                                                      !
!     The Gaussian-distributed random number y is obtained from the    !
!     uniformally ditributed random number x by the following          !
!     equation :                                                       !
!                  y                                                   !
!                 /                                                    !
!          1      [                                                    !
!      ---------- |  exp(- 0.5 * t^2) dt  -  x  =  0                   ! 
!      sqrt(2*pi) ]                                                    !
!                 /                                                    !
!                -inf                                                  !
!                                                                      !
!     We write this equation as F(x,y) = 0. In this program, for       !
!     given x, we solve F(x,y) = 0 by Newton method (we adopt Newton   !
!     method since F(x,y) is monotonous function with respect to y).   !
!     For numerical integration, we adopt Gauss-Legendre formula.      !
!======================================================================!
      module Input
      implicit none
      integer, parameter :: trial_num=10000000
      real*8, parameter :: xmin=-5.0d0, xmax=5.0d0, dx=0.05d0
      end module
!======================================================================!
      program Main
!----------------------------------------------------------------------!
!     In this main program, we make histgram of Gaussian distribution. !
!----------------------------------------------------------------------!
      use Input, only : trial_num, xmin, xmax, dx
      implicit none
      integer :: n, i, j, Ngrid
      real*8, allocatable :: Ny(:)
      real*8, external :: Gaus
      real*8 :: y, x, er
      character(len=15), parameter :: F1='hist.dat'

      Ngrid = anint((xmax - xmin) / dx)
      allocate(Ny(Ngrid))
      Ny = 0.0d0

      do i=1, trial_num
        call Gaus_ran_num(y)
        if (y >= xmin .and. y <= xmax) then
          j = int((y - xmin) / dx) + 1
          Ny(j) = Ny(j) + 1
        end if
      end do

      Ny = Ny / (dble(trial_num) * dx)         ! Normalization
      open(7, file=trim(F1), action='write')

      x = xmin + 0.5d0 * dx
      er = 0.0d0
      do i=1, Ngrid
        write(7,*) x, Ny(i)
        er = er + (Gaus(x) - Ny(i)) / Gaus(x)
        x = x + dx
      end do

      write(6,*) 'Error =', er / dble(Ngrid)
      close(7)

      call Gnu(7, 'gnu', trim(F1))

      stop
      end program
!====================================================================!
      subroutine  Gaus_ran_num(w)
!--------------------------------------------------------------------!
!     This subroutine returns Gaussian-distributed random number w.  !
!--------------------------------------------------------------------!
      implicit none
      real*8, intent(out) :: w
      real*8 :: r

      call unifm_ran_num(r)
      call Newton(r, w)

      return
      end subroutine
!====================================================================!
      subroutine Newton(r, x)
      use Input, only : xmin, xmax, dx
      implicit none
      integer :: i
      integer :: MAXi=200
      real*8, intent(out) :: x
      real*8, intent(in) :: r
      real*8, external :: intgrl
      real*8, parameter :: PI=3.141592653589793d0
      real*8, parameter :: epsa=1.0d-300
      real*8, parameter :: epsr=1.0d-10
      real*8 :: x0
      
      x0 = 0.5d0
      do i=1, MAXi
        x = x0 - intgrl(r,x0) * sqrt(2.0d0*PI) * exp(0.50d0*x0*x0)
        if (abs(x-x0) < epsa + epsr*abs(x)) return
        x0 = x
      end do

      if (i == MAXi) stop 'newton'

      return
      end subroutine
!====================================================================!
      function intgrl(r, w)  result(f)
!--------------------------------------------------------------------!
!     Definition of the function f from which we wil solve an        !
!     equation f(w) = 0.                                             !
!--------------------------------------------------------------------!
      implicit none
      integer :: n
      real*8, intent(in) :: r, w
      real*8 :: f, S, D

      call Integration(n, 0.0d0, w, S, D)
      f = 0.5d0 + S - r

      return
      end function
!====================================================================!
      subroutine unifm_ran_num(x)
      implicit none
      real*8 :: x

      call random_number(x)

      return
      end subroutine
!====================================================================!
      subroutine  Integration(n, a, b, S, D)
!--------------------------------------------------------------------!
!     In this subroutine, we calculate the integral until it         !
!     converges. If 'n' becomes larger than 'MAXn', we display an    !
!     error message and break.                                       !
!--------------------------------------------------------------------!
      use Gauss_Legendre, only : GausLegen
      implicit none
      integer, intent(out) :: n
      integer, parameter :: MAXn=100
      real*8, intent(in) :: a, b
      real*8, intent(out) :: S, D
      real*8, parameter :: epsr=1.0d-10
      real*8, external :: f
      real*8 :: S0

      !S = (b - a) * f(0.5d0*(a+b))   ! Value for n=1.

      n=10
      D=0.0d0
      !do n = 2, MAXn
      !   S0 = S 
         call  GausLegen(n, a, b, S)
      !   D = abs(S - S0)
      !   if (D < epsr * abs(S)) return
      !end do

      !if (n == MAXn) then
      !   write(6,*) 'Calculation did not converge.'
      !   stop
      !end if
      
      return 
      end subroutine
!======================================================================!
      function Gaus(x) result(f)
      implicit none
      real*8, intent(in) :: x
      real*8, parameter :: PI=3.141592653589793d0
      real*8 :: f

      f = exp(- 0.5d0 * x * x) / sqrt(2.0d0 * PI)

      return
      end function
!======================================================================!
      subroutine Gnu(N, Fname, F)
      use Input, only : trial_num
      implicit none
      integer, intent(in) :: N
      character(len=*), intent(in) :: Fname
      character(len=*), intent(in) :: F
      character(len=30) :: C

      write(C,*) trial_num
      open(N, file=Fname)
      write(N,*) 'set term postscript eps enhanced color'
      write(N,*) 'set output "Gaus_ran_num2.eps"'
      write(N,*) 'set size 0.8, 0.8'
      write(N,*) 'set title "Gaussian distribution"'
      write(N,*) 'set xrange [-5:5]'
      write(N,*) 'set yrange [0:0.5]'
      write(N,*) 'set xlabel "X"'
      write(N,*) 'set ylabel "Y"'
      write(N,*) 'set label "Trial : '//trim(adjustl(C))//' times" \'
      write(N,*) 'at -4.5, 0.425'
      write(N,*) 'set xtics 1'
      write(N,*) 'set ytics 0.05'
      write(N,*) 'set grid'
      write(N,*) 'f(x) = exp(-0.5*x*x) / sqrt(2*pi)'
      write(N,*) 'pl "'//F//'"title "Random Number" w boxes, \'
      write(N,*) 'f(x) title  "exp(-0.5*x^2) / sqrt(2{/Symbol p})"'
      write(N,*) 'set term x11'
      close(N)

      return
      end subroutine
