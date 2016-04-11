!======================================================================!
      module global_const
      implicit none
      integer, parameter :: Jmesh = 100
      real(8), parameter :: dt = 0.001d0
      real(8), parameter :: dx = 1.0d0/dble(Jmesh)
      real(8), parameter :: a = dt/(dx*dx)
      end module
!======================================================================!
      program main
      use global_const, only : Jmesh, dt, dx, a
      implicit none
      integer :: i, k
      real(8) :: u(0:Jmesh), Qt(Jmesh-1)
      real(8) :: t, tp, x, P, q, R
      real(8), external :: Exact_sol
      character(len=30) :: c
      character(len=30) :: FM = '(1x,f11.8,2f14.8)'

      q = 1.0d0 + 2.0d0 * a
      P = - a
      R = - a
      Qt(1) = q
      do i=2, Jmesh-1
        Qt(i) = q - P * R / Qt(i-1)
      end do

      call Init(Jmesh, dx, u)
      t = 0.0d0
      tp = 0.0d0
      i = 0
      do while(t < 0.5d0)
        if (t >= tp) then
          write(c,*) i
          c = 'output/'//trim(adjustl(c))//'.dat'
          open(7,file=c)
          write(7,*) '# time =', t
          x = 0.0d0
          do k=0, Jmesh
            write(7,FM) x, u(k), Exact_sol(t,x)
            x = x + dx
          end do
          close(7)
          tp = tp + 0.02d0
          i = i + 1
        end if
        call solve(a, Jmesh, Qt,u)
        t = t + dt
      end do

      stop
      end program
!======================================================================!
      subroutine solve(a, Jmesh, Qt, u)
      implicit none
      integer, intent(in) :: Jmesh
      integer :: i
      real(8), intent(in) :: a, Qt(Jmesh-1)
      real(8), intent(inout) :: u(0:Jmesh)
      real(8) :: P, R, u0(0:Jmesh)

      P = - a
      R = - a

      u0(1) = u(1)
      do i=2, Jmesh-1
        u0(i) = u(i) - P * u0(i-1) / Qt(i-1)
      end do

      u(Jmesh-1) = u0(Jmesh-1) / Qt(Jmesh-1)
      do i=Jmesh-2, 1, -1
        u(i) = (u0(i) - R * u(i+1)) / Qt(i)
      end do

      return
      end subroutine
!======================================================================!
      subroutine Init(Jmesh, dx, u)
      implicit none
      integer, intent(in) :: Jmesh
      integer :: i
      real(8), intent(in) :: dx
      real(8), intent(out) :: u(0:Jmesh)
      real(8), external :: Phi
      real(8) :: x

      x = 0.0d0
      do i=0, Jmesh
        u(i) = Phi(x)
        x = x + dx
      end do

      return
      end subroutine
!======================================================================!
      function Phi(x) result(f)
      implicit none
      real(8), intent(in) :: x
      real(8), parameter :: PI=3.141592653589793d0
      real(8) :: f

      f = sin(PI * x)

      return
      end function
!======================================================================!
      function Exact_Sol(t, x) result(f)
      implicit none
      real(8), intent(in) :: t, x
      real(8), parameter :: PI=3.141592653589793d0
      real(8) :: f

      f = exp(-t*PI*PI) * sin(PI*x)

      return
      end function
!======================================================================!
