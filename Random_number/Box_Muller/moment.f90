!======================================================================!
!     Title  : Box_Muller.f90                                          !
!     Date   : 2009-2-4-Wed                                            !
!     Author : Yusa Shusaku                                            !
!     Last modified : 2009-2-9-Mon                                     !
!                                                                      !
!     A program which generates standard normally distributed random   !
!     numbers from uniformally distributed random numbers by Box-      !
!     Muller method.                                                   !
!======================================================================!
      module Input
      implicit none
      integer, parameter :: trial_num=10000000
      real(8), parameter :: xmin=-5.0d0, xmax=5.0d0, dx=0.05d0
      end module
!======================================================================!
      program main
      use Input
      implicit none
      integer :: i, j
      integer :: idum=-123456789
      real(8) :: y, moment(2), aaa, z
      character(len=15), parameter :: F='hist.dat'

      moment = 0.0d0
      aaa = 0.0d0
      do i=1, trial_num
!       call unifm_ran_num(idum, y)
        call gaus_ran_num(idum, y)
        call gaus_ran_num(idum, z)
        moment(1) = moment(1) + y
        moment(2) = moment(2) + y * y
        aaa = aaa + y * z
      end do
      moment = moment / dble(trial_num)
      aaa = aaa / dble(trial_num)
      write(6,*) '1st moment =', moment(1)
      write(6,*) '2nd moment =', moment(2)
      write(6,*) 'aaa =', aaa

      stop
      end program 
!======================================================================!
      subroutine gaus_ran_num(idum, y)
!----------------------------------------------------------------------!
!     This subroutine transforms uniformally distributed random        !
!     numbers x1 and x2 to standard normally distributed random        !
!     numbers y and ystore by Box-Muller method. We stock the ystore   !
!     until next call of the subroutine.                               !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(inout) :: idum
      real(8), intent(out) :: y
      real(8), parameter :: PI=3.141592653589793d0
      real(8), save :: ystore
      real(8) :: x1, x2, w
      logical, save :: gaus_store=.false.

      if (gaus_store) then
         y = ystore
         gaus_store = .false.
      else
         call unifm_ran_num(idum, x1)
         call unifm_ran_num(idum, x2)
         w = sqrt(-2.0d0 * log(x1))
         y = w * cos(2.0d0 * PI * x2)
         ystore = w * sin(2.0d0 * PI * x2)
         gaus_store = .true.
      end if

      return
      end subroutine
!======================================================================!
      subroutine unifm_ran_num(idum, r)
!----------------------------------------------------------------------!
!     A subroutine which generates a random number uniformally         !
!     distributed in the range [0,1].                                  !
!     The code is taken from "Numerical Recipes in Fortran 90".        !
!----------------------------------------------------------------------!
!     A subroutine which generates a random number uniformally         !
      implicit none
      integer, intent(inout) :: idum
      integer, parameter :: IA=16807,IM=2147483647,IQ=12773,IR=2836
      integer, save :: ix=-1, iy=-1, k
      real(8), intent(out) :: r
      real(8), save :: am

      if (idum <= 0 .or. iy < 0) then
         am = nearest(1.0d0, -1.0d0) / IM
         iy = ior(ieor(888889999, abs(idum)), 1)
         ix = ieor(777755555, abs(idum))
         idum = abs(idum) + 1
      end if
      ix = ieor(ix, ishft(ix,13))
      ix = ieor(ix, ishft(ix,-17))
      ix = ieor(ix, ishft(ix,5))
      k = iy / IQ
      iy = IA * (iy - k * IQ) - IR * k
      if (iy < 0) iy = iy + IM
      r = am * ior(iand(IM, ieor(ix,iy)), 1)

      end subroutine
!======================================================================!
      function  Gaussian(x) result(f)
      implicit none
      real(8), intent(in) :: x
      real(8), parameter :: PI=3.141592653589793d0
      real(8) :: f

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
      write(N,*) 'set output "Box_Muller.eps"'
      write(N,*) 'set size 0.8, 0.8'
      write(N,*) 'set title "Box Muller Method"'
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
      write(N,*) 'pl "'//F//'"title "x1 -> y1" w boxes, \'
      write(N,*) 'f(x) title "exp(-0.5*x^2) / sqrt(2{/Symbol p})"'
      write(N,*) 'set term x11'
      close(N)

      return
      end subroutine
