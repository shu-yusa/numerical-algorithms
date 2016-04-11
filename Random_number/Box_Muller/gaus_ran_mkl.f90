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
      use mkl_vsl_type
      use mkl_vsl
      use Input
      implicit none
      integer :: i, j, Ngrid
      integer :: method, brng, stat
      integer :: idum=123456789
      real(8), allocatable :: Ny(:)
      real(8), external :: Gaussian
      real(8) :: y, x, er, r(trial_num)
      character(len=15), parameter :: F='hist.dat'
      type(vsl_stream_state) :: stream

      Ngrid = nint((xmax - xmin) / dx)
      allocate(Ny(Ngrid))
      Ny = 0.0d0

      brng = vsl_brng_mcg31
      method = vsl_method_dgaussian_boxmuller

      stat = vslnewstream(stream,brng,1) 
      stat = vdrnggaussian(method,stream,trial_num,r,0.0d0,1.0d0)
      stat = vsldeletestream(stream)

      do i=1, trial_num
        if (r(i) >= xmin .and. r(i) <= xmax) then
          j = int((r(i) - xmin) / dx) + 1
          Ny(j) = Ny(j) + 1
        end if
      end do

      Ny = Ny / (dble(trial_num) * dx)        ! Normalization
      open(7, file=trim(F), action='write')
      
      x = xmin + 0.5d0 * dx
      er = 0.0d0
      do i=1, Ngrid
        write(7,*) x, Ny(i)
        er = er + (Gaussian(x) - Ny(i)) / Gaussian(x)
        x = x + dx
      end do
      
      write(6,*) 'Error =', er / dble(Ngrid)
      close(7)
      close(8)

      call Gnu(7, 'gnu', trim(F))
      call system('gnuplot < gnu')

      stop
      end program 
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
