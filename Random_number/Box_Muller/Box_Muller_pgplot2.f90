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
      integer :: pgopen, ios
      integer :: i, j, Ngrid
      integer :: idum = -123456789
      real(4), external :: Gauss
      real(4) :: ymax
      real(4), allocatable :: xr(:)
      real(8) :: y
      real(8), allocatable :: Ny(:)
      character(len=15), parameter :: F='hist.dat'
      character(len=40) :: c

      Ngrid = anint((xmax - xmin) / dx)
      allocate(Ny(Ngrid), xr(Ngrid))
      ymax = 0.5d0 * real(trial_num) * dx

      write(c,*) trial_num
      c = 'Trial : '//trim(adjustl(c))//' times'

      Ny = 0.0d0
      do i=1, trial_num
        call gaus_ran_num(idum, y)
        if (y >= xmin .and. y <= xmax) then
          j = int((y - xmin) / dx) + 1
          Ny(j) = Ny(j) + 1
        end if
      end do

      do i=1, Ngrid
        xr(i) = xmin + (0.5 + real(i-1)) * dx
      end do
      
      ios = pgopen('/xserv')
      call pgsci(9)
      call pgenv(real(xmin), real(xmax), 0.0, ymax,0,1)
      call pgsci(1)
      call pglab('X','Y','Box Muller')
      call pgsci(1)
      call pgtext(-4.5, real(ymax)*0.8, trim(c))
      call pgsci(7)
      call pgbin(Ngrid, xr, real(Ny), .true.)
      call pgsci(8)
      call pgfunx(gauss, Ngrid, real(xmin),real(xmax),1)

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
      function  Gauss(x) result(f)
      use Input, only : trial_num, dx
      implicit none
      real(4), intent(in) :: x
      real(8), parameter :: PI=3.141592653589793d0
      real(4) :: f

      f = exp(- 0.5d0 * x * x) / sqrt(2.0d0 * PI)
      f = real(trial_num) * dx * f

      return
      end function
!======================================================================!
