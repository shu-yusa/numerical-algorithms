!======================================================================!
      module special_fct
      implicit none
      public :: gammaln, arg_gamma, pl, nh, nl, jl
      interface gammaln
        module procedure gammaln_d, gammaln_c
      end interface

      contains
!**********************************************************************!
      function gammaln_d(x) result(g)
      implicit none
      real(8), intent(in) :: x
      real(8), dimension(6) :: coef = (/76.18009172947146d0, &
                  -86.50532032941677d0, 24.01409824083091d0, &
                  -1.231739572450155d0, 0.1208650973866179d-2, &
                  -0.5395239384953d-5/)
      real(8), parameter :: stp = 2.5066282746310005d0
      real(8) :: tmp, g

      if (x < 0) stop 'negative argument in gammaln'
      tmp = x + 5.5d0
      tmp = (x + 0.5d0) * log(tmp) - tmp
      g = tmp + log(stp/x*(1.000000000190015d0 &
         + coef(1)/(x+1.0d0) + coef(2)/(x+2.0d0) + coef(3)/(x+3.0d0) &
         + coef(4)/(x+4.0d0) + coef(5)/(x+5.0d0) + coef(6)/(x+6.0d0)))

      end function
!**********************************************************************!
      function gammaln_c(z) result(g)
      implicit none
      complex(8), intent(in) :: z
      complex(8) :: tmp, g
      real(8), dimension(6) :: coef = (/76.18009172947146d0, &
                  -86.50532032941677d0, 24.01409824083091d0, &
                  -1.231739572450155d0, 0.1208650973866179d-2, &
                  -0.5395239384953d-5/)
      real(8), parameter :: stp = 2.5066282746310005d0

      if (dble(z) < 0) stop 'negative argument in gammaln'
      tmp = z + 5.5d0
      tmp = (z + 0.5d0) * log(tmp) - tmp
      g = tmp + log(stp/z*(1.000000000190015d0 &
         + coef(1)/(z+1.0d0) + coef(2)/(z+2.0d0) + coef(3)/(z+3.0d0) &
         + coef(4)/(z+4.0d0) + coef(5)/(z+5.0d0) + coef(6)/(z+6.0d0)))

      end function
!**********************************************************************!
      function arg_gamma(z) result(t)
      implicit none
      real(8) :: t
      complex(8), intent(in) :: z
      complex(8) :: gam

      gam = exp(gammaln_c(z))
      t = atan(aimag(gam)/dble(gam))
!     t = aimag(gammaln_c(z))
      
      end function
!**********************************************************************!
      pure function PL(L, x) result(P)
!----------------------------------------------------------------------!
!     Legendre polynomials.                                            !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: L
      integer :: i
      real(8), intent(in) :: x
      real(8) :: P0, P1, P

      if (L == 0) then
        P = 1.0d0
        return
      else if (L == 1) then
        P = x
        return
      end if

      P0 = 1.0d0
      P1 = x
      do i=2, L
        P = (dble(2 * i - 1) * x * P1 - dble(i - 1) * P0) / dble(i)
        P0 = P1
        P1 = P
      end do

      return
      end function
!**********************************************************************!
      pure function PL_v(L, x) result(P)
!----------------------------------------------------------------------!
!     Legendre polynomials. Vector version.                            !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: L
      integer :: i
      real(8), intent(in) :: x
      real(8) :: P(0:L)

      P(0) = 1.0d0
      if (L == 0) return
      P(1) = x
      do i=2, L
        P(i) = (dble(2*i-1) * x * P(i-1) - dble(i-1) * P(i-2)) / dble(i)
      end do

      return
      end function
!**********************************************************************!
      pure function NH(n, x) result(f)
!----------------------------------------------------------------------!
!     modified Hermite polynomials.                                    !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: n
      integer :: i
      real(kind=8), intent(in) :: x
      real(kind=8), parameter :: c = 0.7511255444649425d0 ! 1 / PI^(0.25)
      real(kind=8) :: f, f0, f1

      if (n == 0) then
        f = c
        return
      else if (n == 1) then
        f = sqrt(2.0d0) * c * x
        return
      end if

      f0 = c
      f1 = sqrt(2.0d0) * c * x

      do i=1, n-1
        f = (x * f1 - sqrt(dble(i)*0.5d0) * f0) * sqrt(2.0d0/dble(i+1))
        f0 = f1
        f1 = f
      end do

      return
      end function
!**********************************************************************!
      function  jL(L, x)  result(f)
!----------------------------------------------------------------------!
!     sphrical Bessel function.                                        !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: L
      integer :: i
      real(8), intent(in) :: x
      real(8) :: f, j0, j1

      if (x == 0.0d0) then
        if (L == 0) then
          f = 1.0d0
        else
          f = 0.0d0
        end if
        return
      end if

      if (L == 0) then
        f = sin(x) / x
        return
      else if (L == 1) then
        f = (sin(x) - x * cos(x)) / (x * x)
        return
      end if

      j0 = sin(x) / x
      j1 = (sin(x) - x * cos(x)) / (x * x)
      do i=2, L
        f = dble(2 * i - 1) * j1 / x - j0
        j0 = j1
        j1 = f
      end do

      return
      end function
!**********************************************************************!
      function  nL(L, x)  result(f)
!----------------------------------------------------------------------!
!     sphrical Neumann function.                                       !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: L
      integer :: i
      real(8), intent(in) :: x
      real(8) :: f, n0, n1

      if (L == 0) then
        f = - cos(x) / x
        return
      else if (L == 1) then
        f = - (cos(x) + x * sin(x)) / (x * x)
        return
      end if

      n0 = - cos(x) / x
      n1 = - (cos(x) + x * sin(x)) / (x * x)
      do i=2, L
        f = dble(2 * i - 1) / x * n1 - n0
        n0 = n1
        n1 = f
      end do

      return
      end function
!**********************************************************************!
      end module
!======================================================================!
      program main
      use special_fct
      implicit none
      integer :: i, n, ngrid
      integer :: ios, pgopen
      real, allocatable :: xa(:), fa(:,:)
      real(8) :: dx=0.01d0
      real(8) :: x
      real(8), parameter :: PI = 3.141592653589793d0
      complex(8) :: z, c

!     write(6,*) 'factrial'
!     read(5,*) x
!     write(6,*) exp(gammaln(x+1.0d0))

! ** Legendre polynomial **
!     ios = pgopen('/xserv')
!     call pgsch(1.0)
!     call pgenv(-1.0, 1.0, -1.0, 1.0, 0, 1)
!     call pgscf(2)
!     call pglab('x', 'Pn(x)', 'Legendre polynomials')
!     call pgslw(5)
!     ngrid = nint(2.0d0/dx)
!     allocate(xa(ngrid+1),fa(0:6,ngrid+1))
!     do i=1, ngrid+1
!       xa(i) = -1.0 + real(i-1) * dx
!       do n=0, 6
!         fa(n,i) = pl(n,dble(xa(i)))
!       end do
!       fa(0:6,i) = pl_v(6,dble(xa(i)))
!     end do
!     do n=0, 6
!       call pgsci(n+1)
!       call pgtext(0.1, 0.95-0.08*real(n) ,'n='//char(48+n))
!       call pgline(ngrid+1, xa, fa(n,:), 0, 1)
!     end do
!     call pgsci(1)
!     call pgend
!     deallocate(xa,fa)

! ** spherical Bessel **
      ios = pgopen('/xserv')
      call pgsch(1.0)
      call pgenv(0.0, 10.0, -1.0, 1.0, 0, 1)
      call pgscf(2)
      call pglab('x', 'Pn(x)', 'Spherical Bessel')
      call pgslw(5)
      ngrid = nint(10.0d0/dx)
      allocate(xa(ngrid+1),fa(0:6,ngrid+1))
      do i=1, ngrid+1
        xa(i) = real(i-1) * dx
        do n=0, 6
          fa(n,i) = nl(n,dble(xa(i)))
        end do
!       fa(0:6,i) = jl(6,dble(xa(i)))
      end do
      write(6,*) fa(6,1), fa(6,2), fa(6,3),fa(6,4),fa(6,5)
      do n=0, 6
        call pgsci(n+1)
        call pgtext(0.1, 0.95-0.08*real(n) ,'n='//char(48+n))
        call pgline(ngrid+1, xa, fa(n,:), 0, 1)
      end do
      call pgsci(1)
      call pgend
      deallocate(xa,fa)



      stop
      end program

