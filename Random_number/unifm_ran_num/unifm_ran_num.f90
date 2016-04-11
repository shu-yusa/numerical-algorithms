!======================================================================!
!     Title  : unifm_ran_num.f90                                       !
!     Taken from Numerical Recipes                                     !
!                                                                      !
!======================================================================!
      module Input
      implicit none
      integer, parameter :: M=10
      integer, parameter :: trial_num=1000000
      real*8, parameter :: xmin=-5.0d0, xmax=5.0d0, dx=0.05d0
      end module
!======================================================================!
      program Main
      use Input, only : trial_num, xmin, xmax, dx
      implicit none
      integer :: i, j, Ngrid
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

      Ny = Ny / (dble(trial_num) * dx)
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
!======================================================================!
      subroutine Gaus_ran_num(w)
      use Input, only : M
      implicit none
      integer :: i
      integer :: idum =-123456789
      real*8, intent(out) :: w
      real*8 :: r

      w = 0.0d0
      do i=1, M
        call ran(idum, r)
        w = w + sqrt(12.0d0/dble(M)) * (r - 0.5d0)
      end do

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
      write(N,*) 'set output "Gaus_ran_num1.eps"'
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
      write(N,*) 'pl "'//F//'"title "x -> y" w boxes, \'
      write(N,*) 'f(x) title  "exp(-0.5*x^2) / sqrt(2{/Symbol p})"'
      write(N,*) 'set term x11'
      close(N)

      return
      end subroutine
!======================================================================!
      subroutine unifm_ran_num(N, x)
      implicit none
      integer, intent(in) :: N
      integer :: i
      real*8, intent(out) :: x(N)

      do i=1, N
         x(N) = 1.0d0
      !  write(6,*) 'x(i) =', x(i)
      end do

      return
      end subroutine
!======================================================================!
      subroutine ran(idum, r)
      implicit none
      integer, intent(inout) :: idum
      integer, parameter :: IA=16807,IM=2147483647,IQ=12773,IR=2836
      integer, save :: ix=-1, iy=-1, k
      real*8, intent(out) :: r
      real*8, save :: am

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
