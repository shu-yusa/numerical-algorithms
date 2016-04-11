!======================================================================!
!     Title  : GaussLegen.f90                                          !
!     Author : Yusa Shusaku                                            !
!     Date   : 2008-2-29-Fri                                           !
!     Last modified  : 2008-5-25-Sun                                   !      
!                                                                      !
!     A program that computes an integral by using Gauss-Legendre      !
!     formula.                                                         !
!                                                                      !
!     ** Structure of this program **                                  !
!                                                                      !
!     module      Com_var                                              !
!     program     GLegenint                                            !
!     subroutine  Integration                                          !
!     subroutine  GausLegen                                            !
!     subroutine  axis_weight                                          !
!     subroutine  Newton                                               !
!     function    P                                                    !
!     function    f                                                    !
!     subroutine  Draw_Line                                            !
!     subroutine  pict                                                 !
!======================================================================!
      module Com_var
      implicit none
      integer :: MAXn
      real(8) :: epsr, epsa
      parameter(MAXn=100)
      parameter(epsr=1.0d-15, epsa=tiny(epsr))
      end module
!======================================================================!
      program Main
!----------------------------------------------------------------------!
!     Main program.                                                    !
!----------------------------------------------------------------------!
      implicit none
      integer :: n
      real(8)  :: S, err
      real(8)  :: a, b
      character(len=10) :: f100
      parameter(f100='(17x,a)')

      write(6,*)
      write(6,f100) '**************************'    
      write(6,f100) '  Gauss-Legendre formula  '
      write(6,f100) '**************************'    

      a = 0.0d0
      b = 1.0d0
      call Integration(n, a, b, S, err)
      call pict(n, S, err)

      end program
!====================================================================!
      subroutine  Integration(n, a, b, S, D)
!--------------------------------------------------------------------!
!     In this subroutine, we calculate the integral until it         !
!     converges. If 'n' becomes larger than 'MAXn', we display an    !
!     error message and break.                                       !
!--------------------------------------------------------------------!
      use Com_var, only : Maxn, epsr
      implicit none
      integer,intent(out) :: n
      real(8), intent(in) :: a, b
      real(8), intent(out) :: S, D
      real(8), external :: f
      real(8) :: S0

      S = (b - a) * f(0.5d0*(a+b))   ! Value for n=1.

      do n = 2, MAXn
         S0 = S 
         call  GausLegen(n, a, b, S)
         D = abs(S - S0)
         if (D < epsr * abs(S)) return
      end do

      if (n == MAXn) then
         write(6,*) 'Calculation did not converge.'
         stop
      end if
      
      return 
      end subroutine
!====================================================================!
      subroutine  GausLegen(n, a, b, S) 
!--------------------------------------------------------------------!
!     This subroutine performs actual calculation for given 'n'.     !
!     In this subroutine, we transform nodes 't(i)' to the original  !
!     variable 'xp'.                                                 !
!     Since the nodes 't(i)' are symmetric about 't=0.0', it is      !
!     sufficient to calculate only positive nodes. Negative ones are !
!     given by '-t(i)'.                                              !
!--------------------------------------------------------------------!
      implicit none
      integer,intent(in) :: n
      integer :: i, m
      real(8), intent(in) :: a, b
      real(8), intent(out) :: S
      real(8), external :: f
      real(8) :: xp, xn, w(1:n), t(1:n), c, d

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
!====================================================================!
      subroutine  axis_weight(n, m, t, w)
!--------------------------------------------------------------------!
!     In this subroutine, we calculate abscissas and weights of      !
!     Gauss-Legendre formula. Since the nodes of Legendre polynomial !
!     are in the range (-1,1), we adopt cosine function as an        !
!     initial value of Newton method by which we seek for the nodes. !
!--------------------------------------------------------------------!
      implicit none
      integer, intent(in)  :: n
      integer, intent(out) :: m 
      integer :: i 
      real(8), intent(out) :: t(1:n), w(1:n)
      real(8), parameter :: PI=3.141592653589793d0
      real(8)  :: arg, dP, P
      character(len=40) :: f100, f200, f300
      parameter(f100 = '(1x,a,i3,a,f19.15,a,i3,a,f19.15," |")')
      parameter(f200 = '(1x,a,i3,3x,a,5x,"|",10x,a,12x,"|")'  )
      parameter(f300 = '(29x,"|",29x,"|")' )

      arg = PI / dble(4 * n + 2)

      m = (n+1) / 2
!     m = int(0.5d0 * n + 0.51d0)       ! The number of nodes which
                                        ! are not negative.
      call Draw_Line                              
      write(6,f200) 'n =',n,'(abs)Abscissas','Weights'    
      write(6,f300) 

      do i=m, 1, -1
         t(i) = cos(arg * dble(4 * i - 1 ))    ! Initial value.
         call Newton(n, t(i))
         call Legendre(n-1, t(i), P, dP)
         w(i) = 2.0d0 * (1.0d0 - t(i) * t(i)) &
              &  / (dble(n) * P) ** 2
         write(6,f100) 't(',m-i+1,') =',t(i),' | w(',m-i+1,') =',w(i)
      end do

      return
      end subroutine axis_weight
!====================================================================!
      subroutine Newton(n, t)
!--------------------------------------------------------------------!
!     We use this subroutine in seeking for the nodes of Legendre    !
!     polynomial.                                                    !
!--------------------------------------------------------------------!
      use Com_var, only : epsa, epsr
      implicit none
      integer, intent(in)    :: n
      real(8),  intent(inout) :: t
      integer :: j, MAXj
      real(8)  :: t0, dP, P
      parameter(MAXj=200)

      do j=1, MAXj
          t0 = t
          call Legendre(n, t0, P, dP)
          t  = t0 - P / dP
          if (abs(t - t0) < epsa + epsr * abs(t)) return
      end do

      if (j == MAXj) then
          write(6,*) 'We could not obtain the node.'
          stop
      end if
      
      return
      end subroutine
!====================================================================!
      subroutine Legendre(n, x, P, dP)
!--------------------------------------------------------------------!
!     Definition(Calculation) of Legendre polynomial of order 'n'.   !
!     We are using recursion relation in the calculation.            !
!     'dP' is a derivative of Legendre polynomial.                   !
!--------------------------------------------------------------------!
      implicit none
      integer,intent(in) :: n
      integer :: i
      real(8),intent(in)  :: x
      real(8),intent(out) :: dP, P
      real(8)  :: P0, P1

      if (n == 0) then
        P = 1.0d0
        dP = 0.0d0
        return
      end if
      
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
!====================================================================!
      real(8) function f(x)
      implicit none
      real(8),intent(in) :: x

      f = 4.0d0 / (1.0d0 + x * x)

      return
      end function
!====================================================================!
      subroutine  Draw_Line

      write(6,*) '-----------------------------&
                 &----------------------------'
      return
      end subroutine
!====================================================================!
      subroutine pict(n, S, err)
!--------------------------------------------------------------------!
!     Draw a picture.                                                !
!--------------------------------------------------------------------!
      integer, intent(in) :: n
      real(8), intent(in) :: S, err
      character :: Q*17, R*15
      parameter(Q = "(1x,a,1pd22.15,a)")
      parameter(R = '(1x,a,i4,26x,a)')

      call  Draw_Line
      write(6,*) 
      write(6,*) '|**************************************************|'
      write(6,*) '|  The function we calculated.                     |'
      write(6,*) '|                                                  |'
      write(6,*) '|      1                                           |'
      write(6,*) '|     /                                            |'
      write(6,*) '|     [         4                                  |'
      write(6,*) '| S = |  dx -----------  =  pi                     |'
      write(6,*) '|     ]      1  +  x^2                             |'
      write(6,*) '|     /                  =  3.141592653598793...   |'
      write(6,*) '|     0                                            |'
      write(6,*) '|                                                  |'
      write(6,*) '****************************************************'
      write(6,*)
      write(6,*)    '|******************************************|' 
      write(6,*)    '|  Result                                  |'
      write(6,R)    '|    n     =', n  ,                      ' |'
      write(6,Q)    '|    S     =', S  ,              '         |'
      write(6,Q)    '|    ERROR =', err,              '         |'
      write(6,*)    '|                                          |'
      write(6,*)    '********************************************'

      return
      end subroutine
