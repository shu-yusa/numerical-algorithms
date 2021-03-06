!====================================================================!
!     Title  : QR.f90                                                !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-5-21-Wed                                         !
!     Last modified : 2008-5-23-Fri                                  !
!                                                                    !
!     A program which seeks eigenvalues of a given matrix by QR      !
!     method.                                                        !
!     At first, we transform the given matrix to Hessenberg matrix   !
!     by Householder transformation. Then, we determine the shift of !
!     origin and perform QR decomposition. After that, we compute    !
!     the inverse product of the matrix Q and R. By, iterating this  !
!     procedure, we get eigenvalues.                                 !
!                                                                    !
!     *** Structure of this program ***                              !
!                                                                    !
!     module      Global_Constant                                    !
!     program     main                                               !
!                                                                    !
!     subroutine  QR                                                 !
!     subroutine  Householder                                        !
!     subroutine  Orig_Shift                                         !
!     subroutine  QR_Decomp                                          !
!     subroutine  Inv_Product                                        !
!     subroutine  Eigval_2x2                                         !
!                                                                    !
!     subroutine  Line                                               !
!     subroutine  Title                                              !
!     subroutine  Show_Matrix                                        !
!     subroutine  Show_Result                                        !
!====================================================================!
      module  Global_Constant
      implicit none
      real(8), parameter :: epsa = 1.0d-100
      real(8), parameter :: epsr = 1.0d-15
      real(8), parameter :: PI = 3.141592653589793d0
      end module
!====================================================================!
      program main
      implicit none
      integer :: i, j, N, Iter
      real(8), allocatable :: a(:,:), T(:,:), P(:,:)
      complex(8), allocatable :: Eig(:)

!     open(10,file='Matrix_HB_4x4')
!     read(10,*) N
!     allocate(a(N,N), T(N,N), Eig(N), P(N,N))
!     do i=1, N
!         read(10,*) (a(i,j), j=1, N)
!     end do
!     close(10)
      N = 4
      allocate(a(N,N), T(N,N), Eig(N), P(N,N))
      a(1,:) = (/ 4.2d7, -3.6d7,  8.1d13,  1.7d13/)
      a(2,:) = (/-7.3d8,  6.4d8,  2.3d13,  1.8d14/)
      a(3,:) = (/ 6.3d7, -5.6d7, -1.1d14, -2.6d13/)
      a(4,:) = (/-1.1d9,  9.7d8, -3.7d13, -2.7d14/)

      T = a

      call  Title
      write(6,*) '*** Matrix ***'
      call  Show_Matrix(N, a)

      call QR(N, a, P, Eig, Iter)

      write(6,*) '*** Hessenberg Matrix ***'
      T = Matmul(Matmul(Transpose(P), T), P)
      call  Show_Matrix(N, T)
      write(6,*) '*** Upper Triangular Matrix ***'
      write(6,*) '| If there are complex eigenvalues, the matrix  |'
      write(6,*) '| is not exactly triangular, but has block-form |'
      write(6,*) '| parts.                                        |'
      call Show_Matrix(N, a)


      write(6,*) ' Iteration :', Iter, '  Times'
      write(6,*)
      call  Show_Result(N, Eig)

      end program
!====================================================================!
      subroutine  QR(N, a, P, Eig, Iter)
      use Global_Constant
      implicit none
      integer, intent(in) :: N
      integer, intent(out) :: Iter
      integer :: i, m, Neig
      real(8), intent(inout) :: a(N,N)
      real(8), intent(out) :: P(N,N)
      real(8) :: shift, Qk(N,N), errmx
      complex(8), intent(out) :: Eig(N)
      complex(8) :: Eig0(2)

      call  Householder(N, a, P)
      
      m = N
      Neig = 0

      Iter = 0
 AA:  do while(m > 0)
         Iter = Iter + 1
         if (m == 1) then
             Neig = Neig + 1;    Eig(Neig) = a(m,m)
             m = 0
             exit
         end if
         if (m == 2) then
             call Eigval_2x2(N, m, a, Eig0)
             Neig = Neig + 1;    Eig(Neig) = Eig0(1)
             Neig = Neig + 1;    Eig(Neig) = Eig0(2)
             m = 0
             exit
         end if

         call Orig_Shift(N, m, a, shift)
         call QR_Decomp(N, m, a, Qk)
         call Inv_Product(N, m, a, Qk, shift)

         errmx = epsa + epsr * min(abs(a(m-1,m-1)), abs(a(m,m)))
         if (abs(a(m,m-1)) < errmx) then
            Neig = Neig + 1;   Eig(Neig) = cmplx(a(m,m))
            m = m - 1
         end if
         if (m == 2) cycle
         if (abs(a(m-1,m-2)) < errmx) then
            call  Eigval_2x2(N, m, a, Eig0)
            Neig = Neig + 1;    Eig(Neig) = Eig0(1)
            Neig = Neig + 1;    Eig(Neig) = Eig0(2)
            m = m - 2
         end if
      end do AA

      return
      end subroutine
!====================================================================!
      subroutine Householder(N, a, PP)
      implicit none
      integer :: i, k, j
      integer, intent(in) :: N
      real(8), intent(inout) :: a(N,N)
      real(8), intent(out)  :: PP(N,N)
      real(8) :: t, s, c2, w(N), p(N), q(N), ps(N), qs(N), pq(N), cpw

      PP = 0.0d0
      do i=1, N
          PP(i,i) = 1.0d0
      end do

 AA:  do k=1, N-2
          t = 0.0d0
          do i=k+2, N
              t = t + a(i,k) * a(i,k)
          end do
          if (t /= 0.0d0) then
             t = t + a(k+1,k) * a(k+1,k) 
             s = Sign(sqrt(t), a(k+1,k))
             c2 = 1.0d0 / (t + a(k+1,k) * s)
             
             w(k+1) = a(k+1,k) + s
             w(k+2:N) = a(k+2:N,k)
       
             p(1:N) = c2 * Matmul(a(1:N,k+1:N), w(k+1:N))
             q(1:k) = p(1:k)
             cpw = 0.5d0 * c2 * Dot_product(p(k+1:N), w(k+1:N))
             q(k+1:N) = p(k+1:N) - cpw * w(k+1:N)
             ps(k:N) = c2 * Matmul(Transpose(a(k+1:N,k:N)), w(k+1:N))
             qs(1:k) = ps(1:k)
             qs(k+1:N) = ps(k+1:N) - cpw * w(k+1:N)
             
             forall(i=1:N, j=k+1:N) a(i,j) = a(i,j) - q(i) * w(j)
             forall(i=k+1:N, j=k:N) a(i,j) = a(i,j) - w(i) * qs(j)

             pq(1:N) = c2 * Matmul(PP(1:N,k+1:N), w(k+1:N))
             forall(i=1:N, j=k+1:N) PP(i,j) = PP(i,j) - pq(i) * w(j)
          end if
      end do AA
      return
      end subroutine
!====================================================================!
      subroutine  Orig_Shift(N, m, a, shift)
      implicit none
      integer, intent(in) :: N, m
      integer :: i
      real(8), intent(out) :: a(N,N), shift
      real(8) :: b, D2, Eig1, Eig2
      
      b = 0.5d0 * (a(m-1,m-1) + a(m,m))
      D2 = (0.5d0 * (a(m-1,m-1) - a(m,m))) ** 2 + a(m-1,m) * a(m,m-1)

      if (D2 >= 0.0d0) then
         Eig1 = b + sign(sqrt(D2), b)
         Eig2 = (a(m-1,m-1) * a(m,m) - a(m-1,m) * a(m,m-1)) / Eig1
         if (abs(Eig1 - a(m,m)) <= abs(Eig2 - a(m,m))) then
            shift = Eig1
         else
            shift = Eig2
         end if
      else  
         shift = b
      end if

      do i=1, m
          a(i,i) = a(i,i) - shift
      end do

      return
      end subroutine
!====================================================================!
      subroutine  QR_Decomp(N, m, a, Qk)
      implicit none
      integer, intent(in) :: N, m
      integer :: i, j
      real(8), intent(inout) :: a(N,N)
      real(8), intent(out) :: Qk(N,N)
      real(8) :: d, s, c, tau, aji, ajpi, qij, qijp

      Qk(1:m,1:m) = 0.0d0
      do i=1, m
          Qk(i,i) = 1.0d0
      end do

      do j=1, m-1
         d = sqrt(a(j,j) * a(j,j) + a(j+1,j) * a(j+1,j))
         s = a(j+1,j) / d
         c = a(j,j) / d
         if (c == -1.0d0) then
            c = 0.0d0
            s = 0.0d0
         end if
         tau = s / (1.0d0 + c)
      
         a(j,j) = d
         a(j+1,j) = 0.0d0
         do i=j+1, m
            aji  = a(j,i)
            ajpi = a(j+1,i)
            a(j,i)   = a(j,i)   + s * (ajpi - tau * aji)
            a(j+1,i) = a(j+1,i) - s * (aji  + tau * ajpi)
         end do

         do i=1, m
            qij  = Qk(i,j)
            qijp = Qk(i,j+1)
            Qk(i,j)   = Qk(i,j)   + s * (qijp - tau * qij)
            Qk(i,j+1) = Qk(i,j+1) - s * (qij  + tau * qijp) 
         end do
      end do

      return
      end subroutine
!====================================================================!
      subroutine  Inv_Product(N, m, a, Qk, shift)
      implicit none
      integer, intent(in) :: N, m
      integer :: i, L, k
      real(8), intent(in) :: shift, Qk(N,N)
      real(8), intent(inout) :: a(N,N)
      real(8) :: w(N)

      do i=1, m
         w(i:m) = a(i,i:m)
         do k=max(1, i-1), m
            L = min(m,k+1)
            a(i,k) = Dot_product(w(i:L), Qk(i:L,k))
         end do
         a(i,i) = a(i,i) + shift
      end do

      return
      end subroutine
!====================================================================!
      subroutine  Eigval_2x2(N, m, a, Eig)
      implicit none
      integer, intent(in) :: N, m
      real(8), intent(in) :: a(N,N)
      complex(8), intent(out) :: Eig(2)
      real(8) :: b, D2

      b  =  0.5d0 * (a(m-1,m-1) + a(m,m))
      D2 = (0.5d0 * (a(m-1,m-1) - a(m,m))) ** 2 + a(m-1,m) * a(m,m-1)

      if (D2 >= 0.0d0) then
         Eig(1) = b + sign(sqrt(D2), b)
         Eig(2) = (a(m-1,m-1) * a(m,m) - a(m-1,m) * a(m,m-1)) / Eig(1)
      else
         Eig(1) = b + sqrt(- D2) * (0.0d0, 1.0d0)
         Eig(2) = conjg(Eig(1))
      end if

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
      write(6,*) '**************************************************'
      write(6,*) '               Eigenvalue Problem                 '       
      write(6,*) ' After transformation of the matrix to Hessenberg '
      write(6,*) ' matrix, we find engenvalues by QR method.        '
      write(6,*) '**************************************************'
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  Show_Matrix(N, a)
      implicit none
      integer, intent(in) :: N
      integer :: i, j
      real(8), intent(in) :: a(N,N)
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
      complex(8), intent(in) :: t(N)
      character(len=30) :: FM='(2x,es9.2,a,es9.2,a)'

      write(6,*) '*** Eigenvalues ***'
      write(6,*) '-------------------------------'
      do i=1, N
          write(6,FM) dble(t(i)),' + (', dimag(t(i)),') * i'
      end do
      write(6,*)

      return
      end subroutine
