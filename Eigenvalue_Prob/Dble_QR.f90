!====================================================================!
!     Title  : Dble_QR.f90                                           !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-5-23-Fri                                         !
!     Last modified : 2008-5-23-Fri                                  !
!                                                                    !
!     A program which seeks eigenvalues of a given matrix by double  !
!     QR method.                                                     !
!     At first, we transform the given matrix to Hessenberg matrix   !
!     by Householder transformation. Then, we perform double QR      !
!     method. We determine the shift of origin and execute           !
!     Householder transformation. Iterating this procedure the       !
!     matrix becomes diagonal(if there are complex eigenvalues, the  !
!     transformed matrix has block-form parts) and we get            !
!     eigenvalues.                                                   !
!                                                                    !
!     *** Structure of this program ***                              !
!                                                                    !
!     module      Global_Constant                                    !
!     program     main                                               !
!                                                                    !
!     subroutine  Dble_QR                                            !
!     subroutine  Hessenberg                                         !
!     subroutine  Orig_Shift                                         !
!     subroutine  Householder_1                                      !
!     subroutine  Householder_j                                      !
!     subroutine  Eigval_2x2                                         !
!                                                                    !
!     subroutine  Line                                               !
!     subroutine  Title                                              !
!     subroutine  Show_Matrix                                        !
!     subroutine  Show_Result                                        !
!====================================================================!
      module  Global_Constant
      implicit none
      real(8), parameter :: epsa = 1.0d-300
      real(8), parameter :: epsr = 1.0d-15
      real(8), parameter :: PI = 3.141592653589793d0
      end module
!====================================================================!
      program main
      implicit none
      integer :: i, j, N, Iter
      real(8), allocatable :: a(:,:), P(:,:)
      complex(8), allocatable :: Eig(:)

!     open(10,file='Matrix_HB_4x4')
!     read(10,*) N
!     allocate(a(N,N), Eig(N), P(N,N))
!     do i=1, N
!         read(10,*) (a(i,j), j=1, N)
!     end do
!     close(10)

      N = 25
      allocate(a(N,N), Eig(N), P(N,N))
      do i=1, N
         do j=1, i
            a(i,j) = dble(N - i + 1)
         end do
         do j=i+1, N
            a(i,j) = dble(N - j + 1)
         end do
      end do

!     call  Title
!     write(6,*) '*** Matrix ***'
!     call  Show_Matrix(N, a)

      call  Dble_QR(N, a, P, Eig, Iter)

!     write(6,*) '*** Upper Triangular Matrix ***'
!     write(6,*) '| If there are complex eigenvalues, the matrix  |'
!     write(6,*) '| is not exactly triangular, but has block-form |'
!     write(6,*) '| parts.                                        |'
!     call  Show_Matrix(N, a)


      write(6,*) ' Iteration :', Iter, '  Times'
      write(6,*)
      call  Show_Result(N, Eig)

      end program
!====================================================================!
      subroutine  Dble_QR(N, a, P, Eig, Iter)
      use Global_Constant
      implicit none
      integer, intent(in)  :: N
      integer, intent(out) :: Iter
      integer :: j, m, Neig
      real(8), intent(inout) :: a(N,N)
      real(8), intent(out)   :: P(N,N)
      real(8) :: errmx, plus, prod
      complex(8), intent(out) :: Eig(N)
      complex(8) :: Eig0(2), pp(2)

      call  Hessenberg(N, a, P)
      
      m = N
      Neig = 0
      Iter = 0
      pp(1) = (0.0d0, 0.0d0)
      pp(2) = (0.0d0, 0.0d0)
      plus  = 0.0d0
      prod  = 0.0d0

      do while(m > 0)
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

         call  Orig_Shift(N, m, a, pp, plus, prod)
         call  Householder_1(N, m, a, P, plus, prod)
         do j=2, m-1
            call  Householder_j(N, m, j, a, P)
         end do

         errmx = epsa + epsr * min(abs(a(m-1,m-1)), abs(a(m,m)))
         if (abs(a(m,m-1)) < errmx) then
            Neig = Neig + 1;   Eig(Neig) = cmplx(a(m,m))
            pp = (0.0d0, 0.0d0)
            plus = 0.0d0
            prod = 0.0d0
            m = m - 1
         end if

         if (m == 2) cycle
         if (abs(a(m-1,m-2)) < errmx) then
            call  Eigval_2x2(N, m, a, Eig0)
            Neig = Neig + 1;    Eig(Neig) = Eig0(1)
            Neig = Neig + 1;    Eig(Neig) = Eig0(2)
            pp = (0.0d0, 0.0d0)
            plus = 0.0d0
            prod = 0.0d0
            m = m - 2
         end if
      end do

      return
      end subroutine
!====================================================================!
      subroutine Hessenberg(N, a, PP)
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

      do k=1, N-2
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
      end do
      return
      end subroutine
!====================================================================!
      subroutine  Orig_Shift(N, m, a, pp, plus, prod)
      implicit none
      integer, intent(in) :: N, m
      complex(8), intent(inout) :: pp(2)
      real(8), intent(in) :: a(N,N)
      real(8), intent(out) :: plus, prod
      complex(8) :: Eig(2), pm1, pm2
      
      pm1 = pp(1)
      pm2 = pp(2)
      call  Eigval_2x2(N, m, a, Eig)
      pp(1) = Eig(1)
      pp(2) = Eig(2)
      
      if (abs(pp(1) - pm1) <= 0.5d0 * abs(pp(1))) then
         if (abs(pp(2) - pm2) <= 0.5d0 * abs(pp(2))) then
            plus = dble(pp(1) + pp(2))
            prod = dble(pp(1) * pp(2))
         else
            plus = 2.0d0 * dble(pp(1))
            prod = dble(pp(1)) * dble(pp(1))
         end if
      else 
         if (abs(pp(2) - pm2) <= 0.5d0 * abs(pp(2))) then
            plus = 2.0d0 * dble(pp(2)) 
            prod = dble(pp(2)) * dble(pp(2))
         end if
      end if

      return
      end subroutine
!====================================================================!
      subroutine  Householder_1(N, m, a, P, plus, prod)
      use Global_Constant, only: epsa
      implicit none
      integer, intent(in) :: N, m
      integer :: i, k, L
      real(8), intent(inout) :: a(N,N), P(N,N)
      real(8) :: b1, b2, b3, r, r2, c2, cpw, q(N), w(3), qs(N)
      real(8), intent(in) :: plus, prod
       
      b1 = a(1,1) * (a(1,1) - plus) + prod + a(1,2) * a(2,1)
      b2 = a(2,1) * (a(1,1) + a(2,2) - plus)
      b3 = a(3,2) * a(2,1)
      r2 = b2 * b2 + b3 * b3
      if (r2 >= 0.0d0 .or. r2 < 0.0d0) then
      else
         stop 'There is a serious error. Reduce the dimension of matrix'
      end if

      if (r2 /= 0.0d0) then
         r2 = r2 + b1 * b1
         r  = sign(sqrt(r2), b1)
         if (r2 + b1 * r < epsa) return
         c2 = 1.0d0 / (r2 + b1 * r)
         w(1) = b1 + r
         w(2) = b2
         w(3) = b3

         L = min(m,4)
         do i=1, L
            k = max(1,i-1)
            q(i) = c2 * Dot_product(a(i,k:3), w(k:3))
         end do
!        q(1:L) = c2 * Matmul(a(1:L,k:3), w(k:3)) 
         do i=1, L
            k = max(1,i-1)
            q(i) = c2 * Dot_product(a(i,k:3), w(k:3))
         end do
         cpw = 0.5d0 * c2 * Dot_product(q(1:3), w(1:3))
         q(1:3) = q(1:3) - cpw * w(1:3)
         do i=1, m
            k = min(i+1,3)
            qs(i) = c2 * Dot_product(a(1:k,i), w(1:k))
         end do
         qs(1:3) = qs(1:3) - cpw * w(1:3)

         forall(i=1:L, k=1:3)  a(i,k) = a(i,k) - q(i) * w(k)
         forall(i=1:3, k=1:m)  a(i,k) = a(i,k) - w(i) * qs(k)
      
         q(1:N) = c2 * Matmul(P(1:N,1:3), w(1:3))
         forall(i=1:N, k=1:3)  P(i,k) = P(i,k) - q(i) * w(k)
      end if

      return
      end subroutine
!====================================================================!
      subroutine  Householder_j(N, m, j, a, P)
      use Global_Constant, only: epsa
      implicit none
      integer, intent(in) :: N, m, j
      integer :: i, L, S, k
      real(8), intent(inout) :: a(N,N), P(N,N)
      real(8) :: q(N), qs(N), r2, c2, r, w(j:j+2), cpw

      L = min(j+2, m)
      S = min(j+3, m)
      r2 = 0.0d0
      do i=j+1, L
         r2 = r2 + a(i,j-1) * a(i,j-1)
      end do

      if (r2 /= 0.0d0) then
         r2 = r2 + a(j,j-1) * a(j,j-1)
         r = sign(sqrt(r2), a(j,j-1))
         if (r2 + a(j,j-1) * r < epsa) return
         c2 = 1.0d0 / (r2 + a(j,j-1) * r)
         w(j) = a(j,j-1) + r   
         w(j+1:L) = a(j+1:L,j-1)
         
         q(1:L) = c2 * Matmul(a(1:L,j:L), w(j:L))
         if (L < m) then
            q(j+3) = c2 * a(j+3,j+2) * w(j+2)
         end if
         cpw = 0.5d0 * c2 * Dot_product(q(j:L), w(j:L))
         q(j:L) = q(j:L) - cpw * w(j:L)
         
         qs(j-1:m) = c2 * Matmul(Transpose(a(j:L,j-1:m)), w(j:L))
         qs(j:L) = qs(j:L) - cpw * w(j:L)

         forall (i=1:S, k=j:L)   a(i,k) = a(i,k) - q(i) * w(k)
         forall (i=j:L, k=j-1:m) a(i,k) = a(i,k) - w(i) * qs(k)
         
         q(1:N) = c2 * Matmul(P(1:N,j:L), w(j:L))
         forall (i=1:N, k=j:L)  P(i,k) = P(i,k) - q(i) * w(k)
      end if

      return
      end subroutine
!====================================================================!
      subroutine  Eigval_2x2(N, m, a, Eig)
      implicit none
      integer, intent(in) :: N, m
      real(8), intent(in) :: a(N,N)
      complex(8), intent(out) :: Eig(2)
      real(8) :: b, D2, det2x2, Eigw

      b  =  0.5d0 * (a(m-1,m-1) + a(m,m))
      D2 = (0.5d0 * (a(m-1,m-1) - a(m,m))) ** 2 + a(m-1,m) * a(m,m-1)
      det2x2 = a(m-1,m-1) * a(m,m) - a(m-1,m) * a(m,m-1)

      if (D2 >= 0.0d0) then
         Eigw = sign(abs(b) + sqrt(D2), b)
         if (Eigw == 0.0d0) stop '0 divide Eigw'
         if (b > 0.0d0) then
            Eig(1) = det2x2 / Eigw
            Eig(2) = Eigw
         else
            Eig(1) = Eigw
            Eig(2) = det2x2 / Eigw
         end if
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
      character(len=30) :: FM='(2x,f10.6,a,f10.6,a)'

      write(6,*) '*** Eigenvalues ***'
      write(6,*) '-------------------------------'
      do i=1, N
          write(6,FM) dble(t(i)),' + (', dimag(t(i)),') * i'
      end do
      write(6,*)

      return
      end subroutine
