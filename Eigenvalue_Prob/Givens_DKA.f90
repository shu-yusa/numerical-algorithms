!====================================================================!
!     Title  : Givens_DKA.f90                                        !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-5-19-Mon                                         !
!                                                                    !
!     A program which finds eigenvalues of a symmetric matrix.       !
!     At first, we tridiagonalizes the matrix by Givens method and   !
!     after that we find engenvalues by DKA method.                  !
!                                                                    !
!     ***  Structure of this program  ***                            !
!                                                                    !
!     module      Global_Constant                                    !
!     program     main                                               !
!     subroutine  Givens                                             !
!                                                                    !
!     subroutine  DKA                                                !
!     subroutine  Aberth                                             !
!                                                                    !
!     subroutine  Line                                               !
!     subroutine  Title                                              !
!     subroutine  B_sort                                             !
!     function    SOL                                                !
!     subroutine  Exact_Solution                                     !
!     subroutine  Show_Matrix                                        !
!     subroutine  Show_Result                                        !
!                                                                    !
!====================================================================!
      module  Global_Constant
      implicit none
      real*8, parameter :: epsa = 1.0d-100
      real*8, parameter :: epsr = 1.0d-15
      real*8, parameter :: PI = 3.141592653589793d0
      end module
!=====================================================================!
      program main
      implicit none
      integer :: i, j, N, M
      real*8, allocatable :: a(:,:), r(:,:), b(:,:)
      real*8, allocatable :: Diag(:), Sub_Diag(:), Eig(:)
      real*8, allocatable :: T(:,:)
      character(len=30) :: FM='(x,a,i4,a)'

      open(10, file='Matrix_Sym_5x5')
      read(10,*) M
      allocate(b(M,M))
      !allocate(a(N,N), r(N,N), T(N,N))
      !allocate(Diag(N), Sub_Diag(N), Eig(N) )
      do i=1, M
           read(10,*) (b(i,j), j=1, M)
      end do
      close(10)
      
      N = 85
      allocate(a(N,N), r(N,N), T(N,N))
      allocate(Diag(N), Sub_Diag(N), Eig(N) )
      do i=1, N
          do j=1, i
              a(i,j) = dble(N - i + 1)
          end do
          do j=i+1, N
              a(i,j) = dble(N - j + 1)
          end do
      end do
      T = a

      call  Givens(N, a, r)
      call  Title
      write(6,*) '*** Matrix(N = 5) ***'
      call  Show_Matrix(M, b)
      !write(6,*) '*** Transformation Matrix R ***'
      !call  Show_Matrix(N, r)
      !write(6,*) '*** Tridiagonalized Matrix ***'
      !T = Matmul(Matmul(Transpose(r), T), r)
      !call  Show_Matrix(N, T)
      call  Exact_Solution
      write(6,FM) 'Calculation below is for N = ', N, '.'
      write(6,*)

      do i=1, N-1
          Diag(i) = a(i,i)
          Sub_Diag(i) = a(i,i+1)
      end do
      Diag(N) = a(N,N)
      
      call  DKA(N, Diag, Sub_Diag, Eig)
      call  Show_Result(N, Eig)

      end program
!====================================================================!
      subroutine Givens(N, a, r)
      implicit none
      integer :: i, j, p, q
      integer, intent(in) :: N
      real*8, intent(inout) :: a(N,N)
      real*8, intent(out) :: r(N,N)
      real*8 :: s, c, tau, d, dpq, h, g, apj, ajp, ajq, aqj
      real*8 :: rjp, rjq
      
      r = 0.0d0
      do i=1, N
          r(i,i) = 1.0d0
      end do

 AA:   do p=2, N-1
 BB:       do q=p+1, N
              d   = sqrt(a(p-1,p) * a(p-1,p) + a(p-1,q) * a(p-1,q))
              s   = - a(p-1,q) / d
              c   =   a(p-1,p) / d
              tau = s / (1.0d0 + c)
              dpq = 0.5d0 * (a(p,p) - a(q,q))
              h   = 2.0d0 * s * (a(p,q) + s * (dpq    - tau * a(p,q)))
              g   = 2.0d0 * s * (dpq    - s * (a(p,q) + tau * dpq   ))
              
              a(p,p) = a(p,p) - h
              a(q,q) = a(q,q) + h
              a(p,q) = a(p,q) + g
              a(p-1,p) = d
              a(p-1,q) = 0.0d0

              do j=p+1, q-1
                  apj = a(p,j)
                  ajq = a(j,q)
                  a(p,j) = apj - s * (ajq + tau * apj)
                  a(j,q) = ajq + s * (apj - tau * ajq)
              end do
              do j=q+1, N
                  apj = a(p,j)
                  aqj = a(q,j)
                  a(p,j) = apj - s * (aqj + tau * apj)
                  a(q,j) = aqj + s * (apj - tau * aqj)
              end do
              do j=1, N
                  rjp = r(j,p)
                  rjq = r(j,q)
                  r(j,p) = rjp - s * (rjq + tau * rjp)
                  r(j,q) = rjq + s * (rjp - tau * rjq)
              end do
          end do  BB
      end do   AA
    
      return
      end subroutine
!=====================================================================!
      subroutine  DKA(N, Diag, Sub_Diag, Re_Eig)
      use Global_Constant
      implicit none
      integer, intent(in) :: N
      integer :: i, k
      real*8, intent(in)  :: Diag(N), Sub_Diag(N)
      real*8 :: error, err0, err1
      real*8, intent(out) :: Re_Eig(N)  
      complex*16 :: Eig(N)
      complex*16 :: Eig0(N), alpha(N), beta(N)
      complex*16 :: d, d0, d1, dp

      call Aberth(N, Diag, Sub_Diag, Eig)
      alpha = cmplx(Diag)
      beta  = cmplx(Sub_Diag)
      do 
          Eig0 = Eig
          error = 0.0d0
          do i=1, N
              d0 = (1.0d0, 0.0d0)
              d1 = Eig0(i) - alpha(N)
              do k=N-1, 1, -1
                  d = (Eig0(i) - alpha(k)) * d1 &
                    &  - beta(k) * beta(k) * d0
                  d0 = d1
                  d1 = d
              end do
              dp = (1.0d0, 0.0d0)
              do k=1, N
                  if(k == i) cycle
                  dp = dp * (Eig0(i) - Eig0(k))
              end do
              Eig(i) = Eig0(i) - d / dp
              err0 = abs(Eig(i) - Eig0(i))
              err1 = abs(Eig(i)) + abs(Eig0(i))
              error = max(error, err0 / (epsa + epsr * err1))
          end do
          if (error < 1.0d0) exit
      end do
      Re_Eig = dble(Eig)

      return
      end subroutine
!=====================================================================!
      subroutine  Aberth(N, Diag, Sub_Diag, Eig)
      use Global_Constant
      implicit none
      integer, intent(in) :: N
      integer :: i, j
      real*8, intent(in) :: Diag(N), Sub_Diag(N)
      real*8 :: R0, R, S
      complex*16, intent(out) :: Eig(N)
      complex*16 :: ci

      R0 = sum(Diag) / dble(N)
      R = abs(Diag(1) - R0) + abs(Sub_Diag(1))
      do i=2, N-1
          R = max(R, abs(Diag(i) - R0) + &
            & abs(Sub_Diag(i-1)) + abs(Sub_Diag(i)))
      end do
      R = max(R, abs(Diag(N) - R0) + abs(Sub_Diag(N-1)))

      ci = 2.0d0 * PI * (0.0d0, 1.0d0) / dble(N)
      do j=1, N
          Eig(j) = cmplx(R0) &
                 & + R * exp(ci * (dble(j) - 0.75d0))
      end do

      return
      end subroutine
!====================================================================!
      subroutine  Line

      write(6,*) '----------------------------'

      return
      end subroutine
!====================================================================!
      subroutine  Title
      implicit none

      write(6,*)
      write(6,*) '*********************************************'
      write(6,*) ' Eigenvalue Problem of Symmetric Matrix.     '
      write(6,*) ' After tridiagonalizing the matrix by Givens '
      write(6,*) ' method, we find engenvalues by DKA method.  '
      write(6,*) '*********************************************'
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  B_Sort(N, t)
      implicit none
      integer, intent(in) :: N
      integer :: i, k, j
      real*8, intent(inout) :: t(N)
      real*8 :: tmax, temp

      tmax = t(1)
      k = 1
      do j=1, N
          tmax = t(j)
          k = j
          do i=j+1, N
              if (t(i) > tmax) then
                  tmax = t(i)
                  k = i
              end if
          end do
      
          if (tmax /= t(j) ) then
              t(k) = t(j)
              t(j) = tmax
          end if
      end do

      return
      end subroutine
!====================================================================!
      function  SOL(N, k)  result(t)
      use Global_Constant
      implicit none
      integer, intent(in) :: N, k
      real*8 :: t

      t = 0.5d0 / (1.0d0 - &
       & cos((2.0d0*dble(k) - 1.0d0) * PI / (2.0d0*dble(n) + 1.0d0)))

      return
      end function
!====================================================================!
      subroutine  Exact_Solution

      write(6,*) '|** Exact Solution ****************************|'
      write(6,*) '|                                              |'
      write(6,*) '|  N : Dimention of the matrix                 |'
      write(6,*) '|                                              |'
      write(6,*) '|                         1                    |'
      write(6,*) '|  t(k) = -----------------------------------  |' 
      write(6,*) '|          2[ 1 - cos(PI*(2k-1) / (2n+1)) ]    |'
      write(6,*) '|                                              |'
      write(6,*) '|  with k=1, 2, ..., N                         |'
      write(6,*) '|                                              |'
      write(6,*) '************************************************'
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  Show_Matrix(N, a)
      implicit none
      integer, intent(in) :: N
      integer :: i, j
      real*8, intent(in) :: a(N,N)
      character(len=30) :: FM

      FM = '(x,a,' // CHAR(48+N) // 'f7.3,a)'
      do i=1, N
          write(6,FM) ' |', (a(i,j), j=1, N), ' |'
      end do
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  Show_Result(N, t)
      implicit none
      integer, intent(in) :: N
      integer :: i
      real*8, intent(in) :: t(N)
      real*8, external :: SOL
      character(len=30) :: FM='(2x,f11.6,a,1pd13.6)'

      call  B_Sort(N, t)
      write(6,*) '*** Results ***'
      write(6,*) ' Eigenvalue |   Error '
      write(6,*) '-------------------------------'
      do i=1, N
          write(6,FM) t(i),'  | ', abs(SOL(N,i) - t(i))
      end do
      write(6,*)

      return
      end subroutine
!====================================================================!
