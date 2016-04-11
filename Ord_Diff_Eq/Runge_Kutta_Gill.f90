!======================================================================!
!     Title  : Runge_Kutta_Gill.f90                                    !
!     Author : Yusa Shusaku                                            !
!     Date   : 2009-2-19-Thu                                           !
!                                                                      !
!======================================================================!
      program main
      implicit none
      integer :: n, i
      integer, parameter :: Nch=2
      real*8 :: x, y(Nch)
      real*8, parameter :: a=0.0d0, b=1.0d0, h=0.01d0

      n = anint((b - a) / h)
      x = 0.0d0
      y = 1.0d0

      do i=1, n
        write(6,*) x, y(1), y(2)
        x = x + h
        call Runge_Kutta_Gill(Nch, x, y, h)
      end do

      stop
      end program
!======================================================================!
!     subroutine Runge_Kutta_Gill(x, y, h)
!----------------------------------------------------------------------!
!     
!----------------------------------------------------------------------!
!     implicit none
!     integer :: i
!     real*8, intent(in) :: x, h
!     real*8, intent(inout) :: y
!     real*8, parameter :: c1 = 1.0d0 - sqrt(0.5d0)
!     real*8, parameter :: c2 = 1.0d0 + sqrt(0.5d0)
!     real*8, save :: cx(1:4) = (/0.0d0, 0.5d0, 0.0d0, 0.5d0/)
!     real*8, save :: ck(1:4) = (/0.5d0, 1.0d0, 1.0d0, 0.5d0/)
!     real*8, save :: cy(1:4) = (/1.0d0, c1, c2, 0.3333333333333333d0/)
!     real*8, save :: cq(1:4) = (/1.0d0, c1, c2, 1.0d0/)
!     real*8, external :: g
!     real*8 :: xw, z, k, q, r

!     xw = x
!     q = 0.0d0
!     do i=1, 4
!       xw = xw + cx(i) * h
!       k = ck(i) * h * g(xw, y)
!       z = y + cy(i) * (k - q)
!       r = z - y
!       q = q + 3.0d0 * r - cq(i) * k
!       y = z
!     end do

!     return
!     end subroutine
!======================================================================!
      function g(x, y) result(f)
      implicit none
      real*8, intent(in) :: x, y
      real*8 :: f

      f = - x * x * y

      return
      end function
!======================================================================!
      subroutine Runge_Kutta_Gill(n, x, y, h)
      implicit none
      integer, intent(in) :: n
      integer :: i, j
      real(8), intent(in) :: x, h
      real(8), parameter :: c1 = 1.0d0 - sqrt(0.5d0)
      real(8), parameter :: c2 = 1.0d0 + sqrt(0.5d0)
      real(8), save :: cx(1:4) = (/0.0d0, 0.5d0, 0.0d0, 0.5d0/)
      real(8), save :: ck(1:4) = (/0.5d0, 1.0d0, 1.0d0, 0.5d0/)
      real(8), save :: cy(1:4) = (/1.0d0, c1, c2, 0.333333333333333d0/)
      real(8), save :: cq(1:4) = (/1.0d0, c1, c2, 1.0d0/)
      real(8) :: hh, xw
      real(8), intent(inout) :: y(n)
      real(8) :: f(n)
      real(8) :: z, k, q(n), r

      xw = x
      q = 0.0d0
      do j=1, 4
        xw = xw + cx(j) * h
        call fct_RKG(n, xw, y, f)
        hh = ck(j) * h
        do i=1, n
          k = hh * f(i)
          z = y(i) + cy(j) * (k - q(i))
          r = z - y(i)
          q(i) = q(i) + 3.0d0 * r - cq(j) * k
          y(i) = z
        end do
      end do

      return
      end subroutine
!======================================================================!
      subroutine fct_RKG(n, x, y, fct)
      implicit none
      integer, intent(in) :: n
      integer :: i
      real(8), intent(in) :: x, y(n)
      real(8), intent(out) :: fct(n)
      real(8), external :: g

      do i=1, n
        fct(i) = g(x, y(i))
      end do

      return
      end subroutine
