!====================================================================!
!     Title  : GausHerm.f90                                          !
!     Author : Yusa Shusaku                                          ! 
!     Date   : 2008-2-14-Thu                                         !
!     Last modified  : 2008-5-25-Sun                                 !
!                                                                    !
!     A program which calculates an integral by using Gauss-Hermite  !
!     formula. In using this formula, we must know division points   !
!     and weights. Since these values are given in numerical table,  !
!     we can perform integration by giving these values in the       !
!     program. However, it is quite inconvenient to write down       !
!     these values from table every time we perform the integration. !
!     Therefore, we calculate not only the value of integral, but    !
!     also the division points and weights in this program.          !
!                                                                    !
!     ** Structure of this program **                                !
!                                                                    !
!     module      Com_var                                            !
!     program     GHint                                              !
!     subroutine  GausHer                                            !
!     subroutine  Pre_Bisec                                          !
!     subroutine  Bisection                                          !
!     subroutine  zeros                                              !
!     subroutine  weight                                             !
!     function    H                                                  !
!     function    f                                                  !
!     subroutine  Draw_Line                                          !
!     subroutine  Draw_Title                                         !
!     subroutine  pict                                               !
!====================================================================!
      module Com_var
!--------------------------------------------------------------------!
!     Definition of parameters common to several program units.      !
!--------------------------------------------------------------------!
      integer :: mm
      real(8)  :: epsr, epsa, PI, SQRTPI
      parameter(mm=12)
      parameter(epsr=1.0d-15, epsa=1.0d-300)
      parameter(PI=3.1415926535897932d0)
      parameter(SQRTPI=sqrt(PI))
      parameter(MAXn=120)
      end module 
!====================================================================!
      program GHint
!--------------------------------------------------------------------!
!     Main program.                                                  !
!     We calculate the integral until it converges.                  !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer :: n
      real(8) :: S, S0, D
      real(8), external :: f
      character(len=20) :: FM1, FM2, FM3, c
      parameter(FM1='(1x,a,2x,i4,15x,a)')
      parameter(FM2='(1x,a,f20.15)')
      parameter(FM3='(1x,a,1pd24.15)')

!     call Draw_Title
      S = f(0.0d0) * SQRTPI
      open(7,file='weight_gausherm.f90')
      write(7,*) '      subroutine weight(n, m, x, w)'
      write(7,*) '      implicit none'
      write(7,*) '      integer, intent(in) :: n, m'
      write(7,*) '      real(8), intent(out) :: x(m), w(m)'
      write(7,*) '      '
      write(7,*) '      if (n > 120) stop "too large n"'
      write(7,*) '      if (m < (n+1) / 2) stop "too small m"'
      write(7,*)
      write(7,*) '      select case(n)'

      do n=2, MAXn
!       call Draw_Line
        S0 = S
      write(c,*) n
      write(7,*) '        case('//trim(adjustl(c))//')'
        call GausHer(n,S)
        D = abs(S - S0)
!       if (D < epsr * (abs(S) + abs(S0))) exit
        if (n == MAXn) then
!         write(*,*) 'Integration did not converge.'
          exit
        endif
      enddo

      write(7,*) '      end select'
      write(7,*)
      write(7,*) '      return'
      write(7,*) '      end subroutine'
      close(7)

!     call Draw_Line
      call pict
      write(6,*) 
      write(6,FM1) 'm     =',mm, 'Used in the integrand.'
      write(6,FM1) 'n     =',n, 'Computed up to this "n".'
      write(6,FM3) 'S     =',S
      write(6,FM3) 'Error =',D
      write(6,*)  

      end  program
!====================================================================!
      subroutine GausHer(n, S)
!--------------------------------------------------------------------!
!     Calculation of integral using Gauss-Hertmite formula.          !
!     Note that the division points are symmetric about the origin   !
!     of abscissa.                                                   !
!     Input  : n                                                     !
!     Output : S                                                     !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer, intent(in) :: n
      integer :: j, mid
      real(8), intent(out) :: S
      real(8), external :: f
      real(8) :: xj, wj, x(MAXn), wmid, w(MAXn)
      character(len=35) :: f100, f200, f300, c
      character(len=40), parameter :: FM1='(11x,a,f19.15,a)'
      character(len=40), parameter :: FM2='(a,1pd22.15)'
      parameter(f200='("|",2x,a,i3,12x," | ",23x,"|")')
      parameter(f300='("| ",5x,a,6x,"| ",6x,a,10x,"|")')
      parameter(f100='("| ",f19.15," | ",1pd22.15," |")')

      mid  = (n + 1) / 2
      call zeros(n, MAXn, x)
      call weight(n, x(mid), w(mid))

!     write(6,f200) 'n =',n
!     write(6,f300) 'Abscissas','Weights'
      write(c,*) 1
      write(7,FM1,advance='no')'x('//trim(adjustl(c))//') =',x(mid),'d0'
      write(7,FM2)' ; w('//trim(adjustl(c))//') =',w(mid)

      if (mod(n,2) == 1) then
        S = w(mid) * f(x(mid))
      else
        S = 0.0d0
      end if
      do j=mid-1, 1, -1
        call weight(n, x(j), w(j))
!       write(6,f100) x(j), w(j)
      write(c,*) mid - j + 1
      write(7,FM1,advance='no')'x('//trim(adjustl(c))//') =',x(j),'d0'
      write(7,FM2)' ; w('//trim(adjustl(c))//') =',w(j)
        S = S + w(j) * (f(x(j)) + f(-x(j)))
      end do
!     write(6,f100) x(mid), w(mid)

      return
      end subroutine
!====================================================================!
      subroutine  Pre_Bisec(n, x, xLEF, xRI)
!--------------------------------------------------------------------!
!     A subroutine which seeks for startig points of bisection       !
!     method.                                                        !
!     The zeros x(i) of Hermite polynomial of order 'n' are such     !
!     that                                                           !
!         0 <= x(i)^2 <= 4*n + 3                                     !
!     Therefore we seek for the starting points in this reagion.     !
!--------------------------------------------------------------------!
      implicit none
      integer,intent(in) :: n
      real(8), intent(inout) :: x
      real(8), intent(out) :: xLEF, xRI 
      real(8),external :: H
      real(8) :: dH, sgn, dx = 0.1d0

      do
        sgn = sign(1.0d0, H(n,x,dH)) * sign(1.0d0, H(n,x+dx,dH))
        if (sgn < 0.0d0) then
          exit
        else 
          x = x + dx
        end if
        if (x > sqrt(dble(4 * n + 3))) then
          write(6,*) 'Error in ''Pre Bisec'''
          stop
        end if
      end do

      xLEF = x
      xRI  = x + dx

      return
      end  subroutine
!====================================================================!
      subroutine  Bisection(n, x)
!--------------------------------------------------------------------!
!     A subroutine that seeks for zero of Hermite polynomial by      !
!     using bisection method.                                        !
!     Startig points of this method are computed from subroutine     !
!     'pre_bisec'.                                                   !
!     By defining 'MINeps', we can perform the computation well.     !
!     I have not figured out why it is so.                           !
!--------------------------------------------------------------------!
      implicit none
      integer, intent(in)  :: n
      integer :: i
      real(8), intent(out) :: x
      real(8), external :: H
      real(8) :: x0, xLEF, xRI, sgn, dH, MINeps

      call Pre_Bisec(n, x, xLEF, xRI)
      x = 0.0d0
      do
        x0 = x
        x  = 0.5d0 * (xLEF + xRI)
        !MINeps = min(epsa, epsr * (abs(x) + abs(x0)))
        !if ( abs(x - x0) < MINeps) then
        if (abs(x - x0) == 0.0d0) then
          exit
        end if
        sgn = sign(1.0d0, H(n,xLEF,dH)) * sign(1.0d0, H(n,x,dH))
        if (sgn > 0.0d0) then
          xLEF = x
        else
          xRI  = x
        end if
      end do
      
      return
      end subroutine
!====================================================================!
      subroutine  zeros(n, maxn, x)
!--------------------------------------------------------------------!
!     Caluculation of zeros of Hermite polynomial by calling         !
!     subroutine 'bisection'. The array of zeros 'x(i)' are already  !
!     arranged after finding all zeros.                              !
!     Since we seek for zeros succssesively, 'y' is the zero of      !
!     Hermite polynomial after each end of do-loop. Therfore, we     !
!     must shift 'y' slightly forward. Otherwise, the bisection      !
!     subroutine will not work well.                                 !
!     Note that we select the procedure depending on whether 'n' is  !
!     even or not.                                                   !
!--------------------------------------------------------------------!
      implicit none
      logical :: even
      integer, intent(in) :: n, maxn
      integer :: i
      real(8), intent(out) :: x(maxn)
      real(8) :: y, dy = 0.1d0

      even = (mod(n,2) == 0)
      y = 0.0d0
      if (even) then
        do i=(n+1)/2, 1, -1
          y = y + dy
          call Bisection(n, y)
          x(i) = y
          x(n-i+1) = - y
        end do
      else 
        do i=(n-1)/2, 1, -1
          y = y + dy
          call Bisection(n, y)
          x(i) = y
          x(n-i+1) = - y
        end do
        x((n+1)/2) = 0.0d0
      end if

      return
      end subroutine
!====================================================================!
      subroutine weight(n, xj, wj)
!--------------------------------------------------------------------!
!     Calculation of a weight factor 'wj' from a zero point 'xj'.    !
!--------------------------------------------------------------------!
      implicit none
      integer, intent(in)  :: n
      integer :: k
      real(8), intent(in)  :: xj
      real(8), intent(out) :: wj
      real(8), parameter :: SQRTPI = sqrt(3.141592653589793d0)
      real(8), external :: H
      real(8) :: w0, dH
      
      w0 = 0.5d0 * SQRTPI
      do k=1,n
        w0 = w0 * 2.0d0 * dble(k)
      enddo
      wj = w0 / (dble(n) * H(n-1,xj,dH)) ** 2

      return
      end subroutine
!====================================================================!
      function H(n, x, dH) result(f)
!--------------------------------------------------------------------!
!     Definition of Hermite polynomial.                              !
!     To obtain the polynomial of order n, we use the recursion      !
!     relation.                                                      !
!     Derivative function of the Hermite polynomial is also          !
!     calculated.                                                    !
!--------------------------------------------------------------------!
      implicit none
      integer,intent(in) :: n
      real(8), intent(in) :: x
      real(8), intent(out) :: dH
      real(8) :: H0, H1, f
      integer :: j

      if (n == 0) then
        f = 1.0d0
        dH = 0.0d0
        return
      end if
      H1 = 1.0d0
      f  = 2.0d0 * x
      
      do j=2, n
        H0 = H1
        H1 = f
        f = 2.0d0 * (x * H1 - (dble(j) - 1.0d0) * H0)
      enddo
       
      dH = 2.0d0 * dble(n) * H1

      return
      end function
!====================================================================!
      function f(x) result(S)
!--------------------------------------------------------------------!
!     Definition of the test function.                               !
!     We can perform the integration of this function exactly.       !
!     Therefore, we can check whether the program works well or not. !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      real(8), intent(in)  :: x
      real(8) :: S
      integer :: i

      S = 1.0d0
      do i=1, mm
        S = S * (2.0d0 * dble(i) - 1.0d0)
      end do
      S = x ** (2.0d0 * dble(mm)) * 2.0d0 ** dble(mm) / S

      return
      end function
!====================================================================!
      subroutine  Draw_Line

      write(6,*) '-----------------------&
                 &-----------------------'

      end subroutine
!====================================================================!
      subroutine  Draw_Title
      implicit none
      character :: Q*10
      parameter(Q='(9x,a)')

      write(6,*) 
      write(6,Q) '*************************'
      write(6,Q) '  Gauss-Hermite formula  '
      write(6,Q) '*************************'

      end subroutine
!====================================================================!
      subroutine pict
!--------------------------------------------------------------------!
!     This subroutine draws a picture...                             !
!--------------------------------------------------------------------!
      use Com_var
      implicit none
      integer :: i
      real(8) :: T
      character :: Q*25
      parameter(Q='(1x,a,1pd22.15,10x,a)')

      T = 1.0d0
      do i=1, mm
          T = T * (2.0d0 * dble(i) - 1.0d0)
      end do

      write(6,*)
      write(6,*) '|**********************************************|'
      write(6,*) '| The integral we calculated.                  |'
      write(6,*) '|                                              |'
      write(6,*) '|      inf                                     |'
      write(6,*) '|      /                                       |'
      write(6,*) '|      ]  exp(-x^2) 2^m x^(2m)                 |'
      write(6,*) '|  S = | --------------------- dx = sqrt(pi)   |'
      write(6,*) '|      ]       (2m - 1)!!                      |'
      write(6,*) '|      /                                       |'
      write(6,*) '|     -inf                                     |'
      write(6,*) '|                                              |'
      write(6,Q) '|    sqrt(pi) =', SQRTPI,                     '|' 
      write(6,*) '|                                              |'
      write(6,*) '************************************************'

      return 
      end subroutine
!====================================================================!
