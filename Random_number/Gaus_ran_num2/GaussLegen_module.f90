!======================================================================!
      module  Gauss_Legendre
!----------------------------------------------------------------------!
!     --- usage ---                                                    !
!     use Gauss_Legendre, only : GausLegen                             !
!     -- Input parameters --                                           !
!     n    : Number of points of Gauss-Legendre formula                !
!     a, b : Integeration interval                                     !
!     S    : Integrated value                                          ! 
!     -- Note --                                                       !
!     Integrand must be given by another external function. If it has  !
!     two or more parameters, we have to modify this module subroutine.!
!----------------------------------------------------------------------!
      contains
!**********************************************************************!
      subroutine  GausLegen(n, a, b, S) 
!----------------------------------------------------------------------!
!     This subroutine performs actual calculation for given 'n'.       !
!     In this subroutine, we transform nodes 't(i)' to the original    !
!     variable 'xp'.                                                   !
!     Since the nodes 't(i)' are symmetric about 't=0.0', it is        !
!     sufficient to calculate only positive nodes. Negative ones are   !
!     given by '-t(i)'.                                                !
!----------------------------------------------------------------------!
      implicit none
      integer,intent(in)    :: n
      integer :: i, m
      real*8, intent(in) :: a, b
      real*8, intent(out) :: S
      real*8, external :: f
      real*8  :: xp, xn, w(1:n), t(1:n), c, d

      call axis_weight(n, m, t, w)
      S = 0.0d0
      c = 0.5d0 * (b - a)
      d = 0.5d0 * (a + b)

      if (mod(n,2) == 0) then                  ! If n is even...
         do i=1, m
             xp =   c * t(i) + d
             xn = - c * t(i) + d
             S  = S + w(i) * (f(xp) + f(xn))
         end do
      else                                      ! If n is odd...
         do i=1, m-1 
             xp =   c * t(i) + d
             xn = - c * t(i) + d
             S  = S + w(i) * (f(xp) + f(xn))
         end do
         S = S + w(m) * f(d)     ! This point corresponds to
      end if                              ! t(i) = 0.
      
      S = S * c

      return
      end subroutine 
!======================================================================!
      subroutine  axis_weight(n, m, t, w)
!----------------------------------------------------------------------!
!     In this subroutine, we calculate abscissas and weights of        !
!     Gauss-Legendre formula. Since the nodes of Legendre polynomial   !
!     are in the range (-1,1), we adopt cosine function as an          !
!     initial value of Newton method by which we seek for the nodes.   !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in)  :: n
      integer, intent(out) :: m 
      integer :: i 
      real*8, intent(out) :: t(1:n), w(1:n)
      real*8, parameter :: PI=3.141592653589793d0
      real*8  :: arg, dP, P
      character :: f100*40, f200*40, f300*20
      parameter(f100 = '(1x,a,i3,a,f19.15,a,i3,a,f19.15," |")')
      parameter(f200 = '(1x,a,i3,3x,a,5x,"|",10x,a,12x,"|")'  )
      parameter(f300 = '(29x,"|",29x,"|")' )

      arg = PI / dble(4 * n + 2)

      m = int(0.5d0 * n + 0.51d0)       ! The number of nodes which
                                        ! are not negative.
      !call Draw_Line                              
      !write(6,f200) 'n =',n,'(abs)Abscissas','Weights'    
      !write(6,f300) 

      do i=m, 1, -1
         t(i) = cos(arg * dble(4 * i - 1 ))    ! Initial value.
         call Newton(n, t(i))
         call Legendre(n-1, t(i), P, dP)
         w(i) = 2.0d0 * (1.0d0 - t(i) * t(i))                          &
     &            / (dble(n) * P) ** 2
     !     write(6,f100) 't(',m-i+1,') =',t(i),' | w(',m-i+1,') =',w(i)
      end do

      return
      end subroutine 
!======================================================================!
      subroutine Newton(n, t)
!----------------------------------------------------------------------!
!     We use this subroutine in seeking for the nodes of Legendre      !
!     polynomial.                                                      !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: n
      integer, parameter :: jmax=200
      integer :: j
      real*8, intent(inout) :: t
      real*8  :: t0, dP, P
      real*8, parameter :: epsr=epsilon(t0)
      real*8, parameter :: epsa=tiny(t0)

      do j=1, jmax
         t0 = t
         call Legendre(n, t0, P, dP)
         t  = t0 - P / dP
         if ( abs(t - t0) < epsa + epsr * abs(t) ) return
      end do

      return
      end subroutine
!======================================================================!
      subroutine  Legendre(n, x, P, dP)
!----------------------------------------------------------------------!
!     Subroutine(Calculation) for Legendre polynomial of order 'n'.    !
!     We are using recursion relation in the calculation.              !
!     'dP' is a derivative of Legendre polynomial.                     !
!----------------------------------------------------------------------!
      implicit none
      integer,intent(in) :: n
      integer :: i
      real*8,intent(in)  :: x
      real*8,intent(out) :: dP, P
      real*8  :: P0, P1

      P1 = 1.0d0
      P  = x

      do i=2, n
         P0 = P1
         P1 = P
         P = (dble((2*i - 1)) * x * P1 - dble((i - 1)) * P0) / dble(i)
      end do 

      dP = dble(n) * (x * P - P1) / (x * x - 1.0d0)

      return
      end subroutine
!**********************************************************************!
      end module
