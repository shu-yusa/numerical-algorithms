!====================================================================!
!     Title  : Lanczos_GL.f90                                        !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-5-18-Sun                                         !
!     Last modified : 2008-5-19-Mon                                  !
!                                                                    !
!     A program which finds eigenvalues of a symmetric matrix.       !
!     At first, we tridiagonalizes the matrix by Givens method and   !
!     after that we find engenvalues by bisection method.            !
!                                                                    !
!     ***  Structure of this program  ***                            !
!                                                                    !
!     module      Global_Constant                                    !
!     program     main                                               !
!     subroutine  Lanczos                                            !
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
      integer :: i, j, N
      real*8, allocatable :: a(:,:), V(:,:), U(:,:)
      real*8, allocatable :: alpha(:), beta(:)
      complex*16, allocatable :: Eig(:)

      open(10, file='Matrix_HB_4x4')
      read(10,*) N
      allocate(a(N,N), V(N,N), U(N,N))
      allocate(alpha(N), beta(N), Eig(N) )
      do i=1, N
           read(10,*) (a(i,j), j=1, N)
      end do
      close(10)

      call  Title
      write(6,*) '*** Matrix ***'
      call  Show_Matrix(N, a)
      call  Lanczos(N, a, alpha, beta, U, V)
      write(6,*) '*** Transformation Matrix U and V ***'
      call  Show_Matrix(N, U)
      call  Show_Matrix(N, V)

      call  DKA(N, alpha, beta, Eig)
      call  Show_Result(N, Eig)

      end program
!====================================================================!
      subroutine  Lanczos(N, A, alpha, beta, u, v)
      implicit none
      integer, intent(in) :: N
      integer :: k
      real*8, intent(in)  :: A(N,N)
      real*8, intent(out) :: alpha(N), beta(N), v(N,N), u(N,N)
      real*8 :: R1, R2, p(N), q(N), u1(N), v1(N)

      u  = 0.0d0
      u(1,1) = 1.0d0
      v  = u
      R2 = Dot_product(u(:,1), v(:,1))
      p  = Matmul(A, u(:,1))
      q  = Matmul(Transpose(A), v(:,1))
      alpha(1) = Dot_product(v(:,1), p) / R2
      do k=1, N-1
          R1 = R2
          u1 = u(:,k)
          v1 = v(:,k)
          u(:,k+1) = p - alpha(k) * u1
          v(:,k+1) = q - alpha(k) * v1
          R2 = Dot_product(u(:,k+1), v(:,k+1))
          beta(k) = R2 / R1
          p  = Matmul(A, u(:,k+1)) - beta(k) * u1
          q  = Matmul(Transpose(A), v(:,k+1)) - beta(k) * v1
          alpha(k+1) = Dot_product(v(:,k+1), p) / R2
      end do

      return
      end subroutine
!=====================================================================!
      subroutine  DKA(N, Diag, Sub_Diag, Eig)
      use Global_Constant
      implicit none
      integer, intent(in) :: N
      integer :: i, k
      real*8, intent(in)  :: Diag(N), Sub_Diag(N)
      complex*16, intent(out) :: Eig(N)
      complex*16 :: Eig0(N), alpha(N), beta(N)
      complex*16 :: d, d0, d1, dp
      real*8 :: error, err0, err1

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
                  d  = (Eig0(i) - alpha(k)) * d1 - beta(k) * d0
                  d0 = d1
                  d1 = d
              end do
              dp = (1.0d0, 0.0d0)
              do k=1, N
                  if(k == i) cycle
                  dp = dp * (Eig0(i) - Eig0(k))
              end do
              Eig(i) = Eig0(i) - d / dp
              err0   = abs(Eig(i) - Eig0(i))
              err1   = abs(Eig(i)) + abs(Eig0(i))
              error  = max(error, err0 / (epsa + epsr * err1))
          end do
          if (error < 1.0d0) exit
      end do

      return
      end subroutine
!=====================================================================!
      subroutine  Aberth(N, alpha, beta, Eig)
      use Global_Constant
      implicit none
      integer, intent(in) :: N
      integer :: i, j
      real*8, intent(in) :: alpha(N), beta(N)
      real*8 :: R0, R, S
      complex*16, intent(out) :: Eig(N)
      complex*16 :: ci

      R0 = sum(alpha) / dble(N)
      R = abs(alpha(1) - R0) + abs(beta(1))
      do i=2, N-1
          R = max(R, abs(alpha(i) - R0) + 1.0d0 + abs(beta(i)))
      end do
      R = max(R, abs(alpha(N) - R0) + 1.0d0)

      ci = 2.0d0 * PI * (0.0d0, 1.0d0) / dble(N)
      do j=1, N
          Eig(j) = cmplx(R0) &
                 & + R * exp(ci * (dble(j) - 0.75d0))
      end do

      return
      end subroutine
!=====================================================================!
      subroutine  Line

      write(6,*) '----------------------------'

      return
      end subroutine
!====================================================================!
      subroutine  Title
      implicit none

      write(6,*)
      write(6,*) '**********************************************'
      write(6,*) ' Eigenvalue Problem of an Asymmetric Matrix.  '
      write(6,*) ' After tridiagonalizing the matrix by Lanczos '
      write(6,*) ' method, we find engenvalues by DKA method.   '
      write(6,*) '**********************************************'
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  B_Sort(N, t)
      implicit none
      integer, intent(in) :: N
      integer :: i, k, j
      complex*16, intent(inout) :: t(N)
      complex*16 :: temp
      real*8 :: tmax, abs_t(N)

      abs_t = abs(t)
      tmax = abs_t(1)
      k = 1
      do j=1, N
          tmax = abs_t(j)
          k = j
          do i=j+1, N
              if (abs_t(i) > tmax) then
                  tmax = abs_t(i)
                  k = i
              end if
          end do
      
          if (tmax /= abs_t(j)) then
              temp = t(k)
              t(k) = t(j)
              t(j) = temp
              abs_t(k) = abs_t(j)
              abs_t(j) = tmax
          end if
      end do

      return
      end subroutine
!====================================================================!
      subroutine  Show_Matrix(N, a)
      implicit none
      integer, intent(in) :: N
      integer :: i, j
      real*8, intent(in) :: a(N,N)
      character(len=30) :: FM

      FM = "(x,a," // CHAR(48+N) // "f9.3,a)"
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
      complex*16, intent(in) :: t(N)
      character(len=30) :: FM='(2x,2f10.6,a)'

      call  B_Sort(N, t)
      write(6,*) '*** Results ***'
      write(6,*) '       Eigenvalues     '
      write(6,*) '--------------------------'
      do i=1, N
          write(6,FM) t(i),'*i '
      end do
      write(6,*)

      return
      end subroutine
