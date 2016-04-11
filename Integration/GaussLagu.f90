!====================================================================!
!     Title  : GaussLagu.f90                                         !
!     Author : Yusa Shusaku                                          ! 
!     Date   : 2008-3-1-Sat                                          !
!     Last modified : 2008-5-25-Sun                                  !
!                                                                    !
!     A program which calculates an integral by using Gauss-Laguerre !
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
!     pragram     GLint                                              !
!     subroutine  GausLagu                                           !
!     subroutine  Pre_Bisec                                          !
!     subroutine  Bisection                                          !
!     subroutine  zeros                                              !
!     subroutien  weight                                             !
!     function    L                                                  !
!     function    f                                                  !
!     subroutine  pict                                               !
!     subroutine  Draw_Line                                          !
!     subroutine  Draw_Title                                         !
!     subroutine  Draw_Result                                        !
!====================================================================!
      module Com_var
!--------------------------------------------------------------------!
!     Definition of parameters common to several program units.      !
!--------------------------------------------------------------------!
        integer :: MAXn
        real(8)  :: epsa, epsr, PI, SQRTPI
        parameter(MAXn=150)
        parameter(epsr=epsilon(epsr), epsa=tiny(epsa))
        parameter(PI=3.1415926535897932d0)
        parameter(SQRTPI = sqrt(PI))
      end module 
!====================================================================!
      program GLint
!--------------------------------------------------------------------!
!     Main program.                                                  !
!     We calculate the integral until it converges.                  !
!--------------------------------------------------------------------!
      use Com_var, only : SQRTPI, MAXn, epsa, epsr
      implicit none
      integer :: m, n
      real(8) :: S, S0, DS
      real(8), external :: f

      call Draw_Title

      m = 200
      S = f(0.0d0,m) * SQRTPI

      do n=2, MAXn
          call Draw_Line
          S0 = S
          call GausLagu(n, m, S)
          DS = abs(S - S0)
          
          if(DS < epsa + epsr * (abs(S) + abs(S0))) exit
          
          if(n == MAXn) then
              write(*,*) 'Integration did not converge.'
              exit
          endif
      enddo

      call pict
      call Draw_Result(m, n, S, DS)

      end  program
!====================================================================!
      subroutine GausLagu(n, m, S)
!--------------------------------------------------------------------!
!     Calculation of integration by using Gauss-Laguerre formula.    !
!     If logical variable 'useDKA' is 'T' and 'n' is smaller than 12,!
!     we use DKA method for 'n < 12' and bisection method for        !
!     'n > 13'.                                                      !
!     If 'useDKA' is 'F', we use only bisection method throughout    !
!     this program.                                                  !
!     Input  : n                                                     !
!     Output : S                                                     !
!--------------------------------------------------------------------!
      use Com_var, only : MAXn
      implicit none
      integer,intent(in) :: n, m
      real(8),intent(out) :: S
      integer :: j
      real(8)  :: xj, wj, x(MAXn)
      real(8), external :: f
      character(len=35) :: f100,f200,f300
      parameter(f100='("| ",f19.15," | ",1pd22.15," |")')
      parameter(f200='("|",3x,a,i2,12x," | ",23x,"|")')
      parameter(f300='("| ",5x,a,6x,"| ",7x,a,9x,"|")')

      write(6,f200) 'n =', n
      write(6,f300) 'Abscissas','Weights'
      
      call zeros(n, x)

      S = 0.0d0
      do j=n, 1, -1
          call weight(n, x(j), wj)
          write(6,f100) x(j), wj
          S = S + wj * f(x(j),m)
      end do

      return
      end subroutine 
!====================================================================!
      subroutine  Pre_Bisec(n, x, xLEF, xRI)
!--------------------------------------------------------------------!
!     A subroutine which seeks for startig points of bisection       !
!     method.                                                        !
!     The zeros of Laguerre polynomial of order 'n' are such that    !
!         0 <= x(i) <= 4*n - 3                                       !
!     Therefore we seek for the starting points in this reagion.     !
!--------------------------------------------------------------------!
      implicit none
      integer,intent(in)    :: n
      real(8), intent(inout) :: x
      real(8), intent(out)   :: xLEF, xRI 
      real(8)  :: dL, dx=0.05d0
      real(8),external :: L

      do
          if (L(n, x, dL) * L(n, x + dx, dL) < 0.0d0) then
              exit
          else 
              x = x + dx
          end if

          if (x > dble(4 * n - 3)) then
              write(6,*) 'ERROR'
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
!     A subroutine that seeks for zeros of Hermite polynomial by     !
!     using bisection method. Startig points of this method are      !
!     computed from subroutine 'pre_bisec'.                          !
!--------------------------------------------------------------------!
      use Com_var, only : epsr, epsa
      implicit none
      integer,intent(in)  :: n
      real(8), intent(out) :: x
      integer :: i
      real(8)  :: x0, xLEF, xRI
      real(8)  :: sgn, dL
      real(8),external :: L

      call Pre_Bisec(n, x, xLEF, xRI)
      
      x = 0.0d0
      do
          x0 = x
          x  = 0.5d0 * (xLEF + xRI)

          !if (abs(x - x0) < epsr * (abs(x) + abs(x0))) exit
          if (abs(x - x0) == 0.0d0) exit
          
          sgn = sign(1.0d0, L(n,xLEF,dL)) * sign(1.0d0, L(n,x,dL))
          if (sgn > 0.0d0) then
              xLEF = x
          else
              xRI  = x
          end if
      end do
      
      return
      end subroutine
!====================================================================!
      subroutine  zeros(n, x)
!--------------------------------------------------------------------!
!     Caluculation of zeros of Laguerre polynomial by calling        !
!     subroutine 'bisection'. The array of zeros 'x(i)' are already  !
!     arranged after finding all zeros.                              !
!     Since we seek for zeros succssesively, 'y' is the zero of      !
!     Laguerre polynomial after each end of do-loop. Therfore, we    !
!     must shift 'y' slightly forward. Otherwise, the bisection      !
!     subroutine will not work well.                                 !
!--------------------------------------------------------------------!
      use Com_var, only : MAXn
      implicit none
      integer,intent(in) :: n
      real(8),intent(out) :: x(MAXn)
      integer :: i
      real(8)  :: y, dy=0.01d0

      y = 0.0d0
      do i=n, 1, -1
          y = y + dy
          call Bisection(n, y)
          x(i) = y
      end do

      return
      end subroutine
!====================================================================!
      subroutine  weight(n, xj, wj)
!--------------------------------------------------------------------!
!     Calculation of a weight factor 'wj' from a zero point 'xj'.    !
!--------------------------------------------------------------------!
      implicit none
      integer,intent(in)  :: n
      real(8), intent(in)  :: xj
      real(8), intent(out) :: wj
      integer :: k
      real(8)  :: dL
      real(8), external :: L
      
      wj = xj / (dble(n) * L(n-1,xj,dL)) ** 2

      return
      end subroutine
!====================================================================!
      real(8) function L(n, x, dL)
!--------------------------------------------------------------------!
!     Definition of Laguerre polynomial.                             !
!     To obtain the polynomial of order n, we use the recursion      !
!     relation.                                                      !
!     Derivative function of the Laguerre polynomial is also         !
!     calculated.                                                    !
!--------------------------------------------------------------------!
      implicit none
      integer,intent(in) :: n
      real(8),intent(in)  :: x
      real(8)  :: L0, L1, dL
      integer :: j

      L1 = 1.0d0
      L  = - x + 1.0d0
      
      do j=2, n
          L0 = L1
          L1 = L
          L = (- (x - 2.0d0 * dble(j) + 1.0d0) * L1 &
            &  - (dble(j) - 1.0d0) * L0) / dble(j)
      enddo
       
      dL = dble(n) * (L - L1) / x

      return
      end  function
!====================================================================!
      function f(x, m)  result(testf)
!--------------------------------------------------------------------!
!     Definition of the test function.                               !
!     We can perform the integration of this function exactly.       !
!     Therefore, we can check whether the program works well or not. !
!--------------------------------------------------------------------!
      implicit none
      integer, intent(inout) :: m
      integer :: i, n
      real(8), intent(in) :: x
      real(8) :: testf
      
      testf = 1.0d0
      do i=1, m
          testf = testf * x / dble(i)
      end do

      return
      end  function
!====================================================================!
      subroutine pict
!--------------------------------------------------------------------!
!     This subroutine draws a picture...                             !
!--------------------------------------------------------------------!
      call Draw_Line
      write(6,*)
      write(6,*) '|*****************************************|'
      write(6,*) '| The function we calculated.             |'
      write(6,*) '|                                         |'
      write(6,*) '|     inf                                 |'
      write(6,*) '|     /                                   |'
      write(6,*) '|     ]             x^m                   |'
      write(6,*) '| S = | dx exp(-x) ------ =  1            |'
      write(6,*) '|     ]              m!                   |'
      write(6,*) '|     /                                   |'
      write(6,*) '|     0                                   |'
      write(6,*) '|                                         |'
      write(6,*) '*******************************************'

      return 
      end  subroutine
!====================================================================!
      subroutine  Draw_Line
      
      write(6,*) '-----------------------&
                 &-----------------------'

      return 
      end subroutine
!====================================================================!
      subroutine  Draw_Title
      implicit none
      character :: fm*10
      parameter(fm='(8x,a)')

      write(6,*) 
      write(6,fm) '**************************'
      write(6,fm) '  Gauss-Laguerre formula  '
      write(6,fm) '**************************'

      return 
      end subroutine
!====================================================================!
      subroutine  Draw_Result(m, n, S, DS)
      implicit none 
      integer, intent(in) :: n, m
      real(8), intent(in) :: S, DS
      character :: f100*20, f50*15, f10*20
      parameter(f10 ='(1x,a,2x,i3)')
      parameter(f50 ='(1x,a,f20.15)')
      parameter(f100='(1x,a,2x,1pd22.15)')

      write(6,*) 
      write(6,f10)  'm     =',m
      write(6,f10)  'n     =',n
      write(6,f50)  'S     =',S
      write(6,f100) 'Error =',DS
      write(6,*)  

      return
      end subroutine  
