!====================================================================!
!     Title  : Lanczos_Sy_Bi.f90                                     !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-5-18-Sun                                         !
!                                                                    !
!     A program which finds eigenvalues of a symmetric matrix.       !
!     At first, we tridiagonalizes the matrix by Lanczos method and  !
!     after that we find engenvalues by bisection method.            !
!                                                                    !
!     ***  Structure of this program  ***                            !
!                                                                    !
!     module      Global_Constant                                    !
!     program     main                                               !
!     subroutine  Lanczos                                            !
!                                                                    !
!     subroutine  Eigenvalue                                         !
!     subroutine  Bisec                                              !
!     subroutine  Set_Limit                                          !
!     subroutine  Num_Ch_Sign                                        !
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
      real(8), parameter :: epsa = 1.0d-100
      real(8), parameter :: epsr = 1.0d-15
      real(8), parameter :: PI = 3.141592653589793d0
      end module
!=====================================================================!
      program main
      use print_array
      implicit none
      integer :: i, j, N
      real(8), allocatable :: a(:,:), V(:,:)
      real(8), allocatable :: alpha(:), beta(:), Eig(:)

!     open(10, file='Matrix_Sym_5x5')
!     read(10,*) N
      N = 5
      allocate(a(N,N), V(N,N))
      allocate(alpha(N), beta(N), Eig(N) )
!     do i=1, N
!          read(10,*) (a(i,j), j=1, N)
!     end do
!     close(10)
      do i=1, N
         do j=1, i
            a(i,j) = dble(N - i + 1)
         end do
         do j=i+1, N
            a(i,j) = dble(N - j + 1)
         end do
      end do

      call  Title
      write(6,*) '*** Matrix ***'
      call  print_mat(a)
      call  Lanczos(N, a, alpha, beta, V)
      write(6,*) '*** Transformation Matrix V ***'
      call  print_mat(V)
      write(6,*) '*** Tridiagonalized Matrix ***'
      a = Matmul(Matmul(Transpose(V), a), V)
      call  print_mat(a)

      call  Eigenvalue(N, alpha, beta, Eig)
      call  Exact_Solution
      call  Show_Result(N, Eig)

      end program
!====================================================================!
      subroutine  Lanczos(N, A, alpha, beta, v)
      implicit none
      integer, intent(in) :: N
      integer :: k
      real(8), intent(in)  :: A(N,N)
      real(8), intent(out) :: alpha(N), beta(N), v(N,N)
      real(8) :: w(N)

      v = 0.0d0
      v(1,1) = 1.0d0
      w = Matmul(A,v(:,1))
      alpha(1) = Dot_product(v(:,1),w)
      do k=1, N-1
        v(:,k+1) = w - alpha(k) * v(:,k)
        beta(k)  = sqrt( Dot_product(v(:,k+1), v(:,k+1)) )
        v(:,k+1) = v(:,k+1) / beta(k)
        w = Matmul(A,v(:,k+1)) - beta(k) * v(:,k)
        alpha(k+1) = Dot_product(v(:,k+1),w)
      end do

      return
      end subroutine
!=====================================================================!
      subroutine  Eigenvalue(N, Diag, Sub_Diag, Eig)
      use Global_Constant
      implicit none
      integer, intent(in) :: N
      integer :: Nsl, Nsr, m
      integer :: i, Nl, Nr, Nm
      real(8) :: a, b, Xsr, x, Xsl
      real(8), intent(out) :: Eig(N)
      real(8), intent(in)  :: Diag(N), Sub_Diag(N)
      real(8) :: Xr, Xl, Sub_Diag2(N), Xm, NCS_a, Xl0

      call  Set_Limit(N, Diag, Sub_Diag, Sub_Diag2, a, b)
      m = 0
      Xsl = a
      Xsr = b
      call  Num_Ch_Sign(N, Diag, Sub_Diag2, a, Nsl)
      call  Num_Ch_Sign(N, Diag, Sub_Diag2, b, Nsr)
      NCS_a = Nsl
       
      do
          Xl = Xsl
          Xr = Xsr
          Nl = Nsl 
          Nr = Nsr
          if (Nsl - Nsr == 1) then
              Xl0 = Xl
              call  Bisec(N, Diag, Sub_Diag2, Xl, Xr, Nl, Nr, x)
              m = m + 1
              EIG(m) = x
              Xr = Xl0
              Xl = a
              Nr = Nl
              Nl = NCS_a
              if (m == N) return
          end if
           
          Xm = 0.5d0 * (Xl + Xr)
          call  Num_Ch_Sign(N, Diag, Sub_Diag2, Xm, Nm)

          if (Nl > Nm) then
              Xsl = Xl
              Xsr = Xm
              Nsl = Nl
              Nsr = Nm
          end if
          if (Nm > Nr) then
              Xsl = Xm
              Xsr = Xr
              Nsl = Nm
              Nsr = Nr
          end if
      end do

      return
      end subroutine
!=====================================================================!
      subroutine  Bisec(N, Diag, Sub_Diag2, Xl, Xr, Nl, Nr, x)
      use Global_Constant
      implicit none
      integer :: Nm, N
      integer, intent(in)    :: Nl, Nr
      real(8),  intent(inout) :: Xl, Xr
      real(8),  intent(in)    :: Diag(N), Sub_Diag2(N)
      real(8),  intent(out)   :: x
      real(8) :: eps
        
      do
          x = 0.5d0 * (Xl + Xr)
          call  Num_Ch_Sign(N, Diag, Sub_Diag2, x, Nm)
          if (Nm == Nr) then
              Xr = x
          else if (Nm == Nl) then
              Xl = x
          end if
          eps = epsa + epsr * (abs(Xr) + abs(Xl))
          if (Xr - Xl < eps) exit
      end do
      
      return
      end subroutine
!=====================================================================!
       subroutine  Set_Limit(N, Diag, Sub_Diag, Sub_Diag2, a, b)
       implicit none
       integer, intent(in) :: N
       integer :: i
       real(8), intent(in) :: Diag(N), Sub_Diag(N)
       real(8), intent(out) :: Sub_Diag2(N), a, b

       a = Diag(1) - abs(Sub_Diag(1))
       b = Diag(1) + abs(Sub_Diag(1))
       
       do i=2, N-1
           a = Min(a, Diag(i) - abs(Sub_Diag(i-1)) - abs(Sub_Diag(i)))
           b = Max(b, Diag(i) + abs(Sub_Diag(i-1)) + abs(Sub_Diag(i)))
           Sub_Diag2(i-1) = Sub_Diag(i-1) * Sub_Diag(i-1)
       end do

       a = Min(a, Diag(N) - abs(Sub_Diag(N-1)))
       b = Max(b, Diag(N) + abs(Sub_Diag(N-1)))
       Sub_Diag2(N-1) = Sub_Diag(N-1) * Sub_Diag(N-1)

       return
       end subroutine
!=====================================================================!
       subroutine  Num_Ch_Sign(N, Diag, Sub_Diag2, x, NCS)
       implicit none
       integer :: N, i, M
       integer, intent(inout) :: NCS
       real(8), intent(in) :: Diag(N), Sub_Diag2(N)
       real(8) :: G, G0
       real(8), intent(in) :: x
       real(8), parameter :: eps = 1.0d-10
       
       G0 = x - Diag(1)

       if (G0 < 0.0d0) then
          NCS = 1      
       else 
          NCS = 0
       end if

       do i=2, N
           if (G0 == 0.0d0) G0 = eps
           G = x - Diag(i) - Sub_Diag2(i-1) / G0
           if (G < 0.0d0) NCS = NCS + 1
           G0 = G
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
      write(6,*) '*************************************************'
      write(6,*) ' Eigenvalue Problem of Symmetric Matrix.         '
      write(6,*) ' After tridiagonalizing the matrix by Lanczos    '
      write(6,*) ' method, we find engenvalues by bisection method '
      write(6,*) '*************************************************'
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  B_Sort(N, t)
      implicit none
      integer, intent(in) :: N
      integer :: i, k, j
      real(8), intent(inout) :: t(N)
      real(8) :: tmax, temp

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
      real(8) :: t

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
      real(8), intent(in) :: t(N)
      real(8), external :: SOL
      character(len=30) :: FM='(2x,f9.6,a,1pd13.6)'

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
