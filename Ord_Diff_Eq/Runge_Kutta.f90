!======================================================================!
      program main
      implicit none
      integer :: t1, t2
      real(8) :: xp, x, k
      real(8), dimension(2) :: y, y0
      real(8), parameter :: h = 0.001d0
      character(len=30), parameter :: FM = '(1x,2f11.5,1pd14.5)'

      call system_clock(t1)
      xp = 0.0d0

      x = 0.0d0
      k = 2.0d0
      call Init(x, y)
      do while (x <= 5000.0d0)
        if (x > xp) then
          write(6,FM) x, y(1), abs(y(1) - cos(k*x))
          xp = xp + 100.0d0
        end if
        call runge_kutta(k, x, y, h)
        x = x + h
      end do
      call system_clock(t2)
      write(6,*) 'time =', t2 - t1

      stop
      end program
!======================================================================!
      subroutine Init(x, y)
      implicit none
      real(8), intent(in) :: x
      real(8), intent(out) :: y(2)

      y(1) = 1.0d0
      y(2) = 0.0d0

      return
      end subroutine
!======================================================================!
      subroutine runge_kutta(k, x, y, h)
      implicit none
      real(8), intent(in) :: k, x, h
      real(8), intent(inout) :: y(2)
      real(8), dimension(2) :: k1, k2, k3, k4
      real(8) :: h2

      h2 = 0.5d0 * h
      call fct(k, x, y, k1)
      call fct(k, x+h2, y+h2*k1, k2)
      call fct(k, x+h2, y+h2*k2, k3)
      call fct(k, x+h,  y+h*k3,  k4)
      y = y + h * (k1 + 2.0d0 * (k2 + k3) + k4) / 6.0d0

      return
      end subroutine
!======================================================================!
      subroutine fct(k, x, y, f)
      implicit none
      real(8), intent(in) :: k, x, y(2)
      real(8), intent(out) :: f(2)

      f(1) = y(2)
      f(2) = - k * k * y(1)

      return
      end subroutine
