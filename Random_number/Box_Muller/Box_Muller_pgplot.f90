!======================================================================!
!     Title  : Box_Muller_pgplot.f90                                   !
!     Date   : 2009-2-4-Wed                                            !
!     Author : Yusa Shusaku                                            !
!     Last modified : 2009-5-17-Sun                                    !
!                                                                      !
!     A program which generates standard normally distributed random   !
!     numbers from uniformally distributed random numbers by Box-      !
!     Muller method.                                                   !
!======================================================================!
      module Input
      implicit none
      integer, parameter :: trial_num=5000000
      real(8), parameter :: xmin=-5.0d0, xmax=5.0d0, dx=0.05d0
      end module
!======================================================================!
      program main
      use Input
      implicit none
      integer :: pgopen, ios
      integer :: i, Ngrid, rgrid
      integer(8) :: idum=-123456789
      real(8), allocatable :: Ny(:)
      real(8), external :: Gauss
      real(8) :: y
      real(4) :: dat(trial_num), ymax
      real(4), allocatable, dimension(:) :: xr, yr
      character(len=40) :: c

      Ngrid = nint((xmax - xmin) / dx)
      rgrid = nint((xmax - xmin) / 0.01d0)
      allocate(xr(rgrid), yr(rgrid))
      ymax = 0.5d0 * dble(trial_num) * dx

      write(c,*) trial_num
      c = 'Trial : '//trim(adjustl(c))//' times'

      do i=1, trial_num
        call gaus_ran_num(idum, y)
        dat(i) = real(y)
      end do

      do i=1, rgrid
        xr(i) = xmin + real(i-1) * 0.01
        yr(i) = Gauss(dble(xr(i))) * ymax * 2.0
      end do

      ios = pgopen('/xserv')
      call pgenv(real(xmin),real(xmax),0.0,ymax,0,1)
      call pglab('X','Y','Box Muller Method')
      call pgtext(-4.5,real(ymax)*0.8,trim(c))
      call pgsci(3)
      call pghist(trial_num, dat, real(xmin), real(xmax), Ngrid, 1)
      call pgsci(7)
      call pgline(rgrid,xr,yr,0,1)

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
      implicit none
      real(8), intent(in) :: x
      real(8), parameter :: PI=3.141592653589793d0
      real(8) :: f

      f = exp(- 0.5d0 * x * x) / sqrt(2.0d0 * PI)

      return
      end function
