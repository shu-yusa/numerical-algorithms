!====================================================================!
!     Title  : Houshold_Bisec.f90                                    !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-5-17-Sat                                         !
!     Last modified : 2009-5-9-Sat                                   !
!                                                                    !
!     A program which finds engenvalues of a symmetric matrix.       !
!     At first, we tridiagonalize the matrix by Householder          !
!     transformation, then we find engenvalues by bisection method.  !
!                                                                    !
!     ***  Structure of this program  ***                            !
!                                                                    !
!     module      Global_Constant                                    !
!     program     main                                               !
!     subroutine  Householder                                        !
!                                                                    !
!     subroutine  Eigenvalue                                         !
!     subroutine  Bisec                                              !
!     subroutine  Set_Limit                                          !
!     subroutine  Num_Ch_Sign                                        !
!     subroutine  Eigenvector                                        !
!     function    Vec_Norm                                           !
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
      module Eigen
!----------------------------------------------------------------------!
!     This program is a module subroutine which seeks eigenvalues      !
!     (and eigenvectors if you want) by Householer-Bisection method.   !
!     ** USAGE **                                                      !
!     use Eigen, only : Houshold_Bisec                                 !
!     -- Input parameters --                                           !
!     N       : Dimension of the matrix                                !
!     A       : Matrix                                                 !
!     vec     : Logical variable. If you want to get eigen vectors,    !
!               vec must be .true..                                    !
!     epsr    : Torrerable relative error                              !
!     Eig_val : Eigen values                                           !
!     Eig_vec : Eigen vectors                                          !
!----------------------------------------------------------------------!
      contains
!**********************************************************************!
      subroutine Household_Bisec(N, A, vec, epsr, Eig_val, Eig_vec)
      implicit none
      integer, intent(in) :: N
      real(8), intent(inout) :: A(N,N)
      real(8), intent(in) :: epsr
      real(8), intent(out) :: Eig_val(N), Eig_vec(N,N)
      real(8), allocatable :: P(:,:)
      real(8), dimension(N) :: alpha, beta
      logical, intent(in) :: vec

      allocate(P(N,N))
      call Householder(N, A, alpha, beta, vec, P)
      call Eigenvalue_up(N, alpha, beta, Eig_val, epsr)
      if (vec) then
        call Eigenvector(N, alpha, beta, Eig_val, P, Eig_vec)
      end if
      deallocate(P)

      return
      end subroutine
!======================================================================!
      subroutine Householder(N, a, alpha, beta, vec, PP)
!----------------------------------------------------------------------!
!     alpha(1:N)  ----- diagonal elements                              !
!     beta(1:N-1) ----- subdiagonal elements                           !
!----------------------------------------------------------------------!
      implicit none
      logical, intent(in) :: vec
      integer :: i, k
      integer, intent(in) :: N
      real(8), intent(inout) :: a(N,N)
      real(8), intent(out) :: alpha(N), beta(N), PP(N,N)
      real(8), allocatable :: PH(:,:)
      real(8), dimension(N) :: p, q, w
      real(8) :: t, s, c2

      if (vec) then
        allocate(PH(N,N))
        PP = 0.0d0
        do i=1,N
          PP(i,i) = 1.0d0
        end do
      end if

      do k=1, N-2
        t = dot_product(a(k+1:N,k),a(k+1:N,k))
        s = Sign(sqrt(t), a(k+1,k))
        alpha(k)  = a(k,k)
        beta(k) = - s

        if (t /= 0.0d0) then
          if (t + a(k+1,k) * s == 0.0d0) stop '0 divide'
            c2 = 1.0d0 / (t + a(k+1,k) * s)
            w(k+1) = a(k+1,k) + s
          else 
            write(6,*) 't = 0'
            cycle
        end if
        w(k+2:N) = a(k+2:N,k)
         
        do i=k+1, N
          p(i) = c2 * (dot_product(a(i,k+1:i), w(k+1:i))     &
     &              +  dot_product(a(i+1:N,i), w(i+1:N)))
        end do

        q(k+1:N) = p(k+1:N) - 0.5d0 * c2 * w(k+1:N)          &
     &            * dot_product(p(k+1:N), w(k+1:N))

        do i=k+1, N
          a(i:N,i) = a(i:N,i) - w(i:N) * q(i) - q(i:N) * w(i)
        end do

        if (vec) then
          PH = 0.0d0
          do i=1, N
            PH(i,i) = 1.0d0
          end do
          do i=k+1, N
            PH(i:N,i) = PH(i:N,i) - c2 * w(i) * w(i:N)
            PH(i,i+1:N) = PH(i+1:N,i)
          end do
          PP = matmul(PP, PH)
        end if
      end do

      alpha(N-1) = a(N-1,N-1)
      alpha(N)   = a(N,N)
      beta(N-1)  = a(N,N-1)
      if (vec) deallocate(PH)

      return
      end subroutine
!======================================================================!
      subroutine  Eigenvalue_down(N, Diag, Sub_Diag, Eig, epsr)
      implicit none
      integer, intent(in) :: N
      integer :: Nsl, Nsr, m
      integer :: Nl, Nr, Nm
      real(8), intent(in)  :: Diag(N), Sub_Diag(N), epsr
      real(8), intent(out) :: Eig(N)
      real(8) :: a, b, Xsr, x, Xsl, G
      real(8) :: Xr, Xl, Sub_Diag2(N), Xm, NCS_a, Xl0

      call  Set_Limit(N, Diag, Sub_Diag, Sub_Diag2, a, b)
      m = 0
      Xsl = a
      Xsr = b
      call  Num_Ch_Sign(N, Diag, Sub_Diag2, a, Nsl, G)
      if (G == 0.0d0) then
        Eig(N) = a
        m = m + 1
      end if
      call  Num_Ch_Sign(N, Diag, Sub_Diag2, b, Nsr, G)
      if (G == 0.0d0) then
        Eig(1) = b
        m = m + 1
      end if
      NCS_a = Nsl
       
      if (m == N) return
      do
        Xl = Xsl ; Xr = Xsr
        Nl = Nsl ; Nr = Nsr
        if (Nsl - Nsr == 1) then
          Xl0 = Xl
          call  Bisec(N, Diag, Sub_Diag2, Xl, Xr, Nl, Nr, x, epsr)
          m = m + 1
          EIG(m) = x
          Xr = Xl0
          Xl = a
          Nr = Nl
          Nl = NCS_a
          if (m == N) return
        end if
         
        Xm = 0.5d0 * (Xl + Xr)
        call  Num_Ch_Sign(N, Diag, Sub_Diag2, Xm, Nm, G)

        if (Nl > Nm) then
          Xsl = Xl ; Xsr = Xm
          Nsl = Nl ; Nsr = Nm
        end if
        if (Nm > Nr) then
          Xsl = Xm ; Xsr = Xr
          Nsl = Nm ; Nsr = Nr
        end if
      end do

      return
      end subroutine
!======================================================================!
      subroutine  Eigenvalue_up(N, Diag, Sub_Diag, Eig, epsr)
      implicit none
      integer, intent(in) :: N
      integer :: Nsl, Nsr, m
      integer :: Nl, Nr, Nm
      real(8), intent(in)  :: Diag(N), Sub_Diag(N), epsr
      real(8), intent(out) :: Eig(N)
      real(8) :: a, b, Xsr, x, Xsl, G
      real(8) :: Xr, Xl, Sub_Diag2(N), Xm, NCS_b, Xr0

      call  Set_Limit(N, Diag, Sub_Diag, Sub_Diag2, a, b)
      m = 0
      Xsl = a
      Xsr = b
      call  Num_Ch_Sign(N, Diag, Sub_Diag2, a, Nsl, G)
      if (G == 0.0d0) then
        Eig(1) = a
        m = m + 1
      end if
      call  Num_Ch_Sign(N, Diag, Sub_Diag2, b, Nsr, G)
      if (G == 0.0d0) then
        Eig(N) = b
        m = m + 1
      end if
      NCS_b = Nsr
       
      if (m == N) return
      do
        Xl = Xsl ; Xr = Xsr
        Nl = Nsl ; Nr = Nsr
        if (Nsl - Nsr == 1) then
          Xr0 = Xr
          call  Bisec(N, Diag, Sub_Diag2, Xl, Xr, Nl, Nr, x, epsr)
          m = m + 1
          EIG(m) = x
          Xr = b
          Xl = Xr0
          Nl = Nr
          Nr = NCS_b
          if (m == N) return
        end if
        Xm = 0.5d0 * (Xl + Xr)
        call  Num_Ch_Sign(N, Diag, Sub_Diag2, Xm, Nm, G)
        if (Nm > Nr) then
          Xsl = Xm ; Xsr = Xr
          Nsl = Nm ; Nsr = Nr
        end if
        if (Nl > Nm) then
          Xsl = Xl ; Xsr = Xm
          Nsl = Nl ; Nsr = Nm
        end if
      end do

      return
      end subroutine
!======================================================================!
      subroutine  Bisec(N, Diag, Sub_Diag2, Xl, Xr, Nl, Nr, x, epsr)
      implicit none
      integer, intent(in)    :: Nl, Nr
      integer :: Nm, N
      real(8), intent(inout) :: Xl, Xr
      real(8), intent(in)    :: Diag(N), Sub_Diag2(N), epsr
      real(8), intent(out)   :: x
      real(8), parameter :: epsa=1.0d-300
      real(8) :: G
        
      do
        x = 0.5d0 * (Xl + Xr)
        call  Num_Ch_Sign(N, Diag, Sub_Diag2, x, Nm, G)
        if (Nm == Nr) then
          Xr = x
        else if (Nm == Nl) then
          Xl = x
        end if
        if (Xr - Xl < epsa + epsr * (abs(Xr) + abs(Xl))) exit
      end do
      
      return
      end subroutine
!======================================================================!
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
!======================================================================!
       subroutine Num_Ch_Sign(N, Diag, Sub_Diag2, x, NCS, G)
       implicit none
       integer :: N, i
       integer, intent(inout) :: NCS
       real(8), intent(in) :: Diag(N), Sub_Diag2(N), x
       real(8), intent(out) :: G
       real(8), parameter :: eps = 1.0d-10
       real(8) :: G0
       
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
!======================================================================!
      subroutine  Eigenvector(N, alpha, beta, Eig, P, x)
      implicit none
      integer, intent(in) :: N
      integer :: i, k
      real(8), intent(in) :: alpha(N), beta(N), Eig(N), P(N,N)
      real(8), intent(out) :: x(N,N)
      real(8) :: v(N)

      do i=1, N
        v(1) = 1.0d0
        if ( beta(1) /= 0.0d0 ) then
          v(2) = - (alpha(1) - Eig(i)) * v(1) / beta(1)
        end if
        do k=2, N-1
          v(k+1) = - (beta(k-1) * v(k-1)                             &
     &                + (alpha(k) - Eig(i)) * v(k)) / beta(k)
        end do
        x(:,i) = Matmul(P, v) / sqrt(dot_product(v,v))
      end do

      return
      end subroutine 
!**********************************************************************!
      end module
!======================================================================!
      program main
      use Eigen, only : Household_Bisec
      implicit none
      integer :: i, j, N
      real(8), allocatable :: a(:,:)
      real(8), allocatable :: Eig(:), x(:,:)
      real(8), parameter :: epsr=1.0d-15
      logical, parameter :: vec=.true.

      N = 200
      allocate(a(N,N))
      allocate(Eig(N), x(N,N))
      do i=1, N
         do j=1, i
            a(i,j) = dble(N - i + 1)
         end do
         do j=i+1, N
            a(i,j) = dble(N - j + 1)
         end do
      end do
!     N = 2
!     allocate(a(2,2))
!     allocate(eig(2), x(2,2))
!     a(1,1) = 0.0d0
!     a(1,2) = 0.25d0
!     a(2,1) = a(1,2)
!     a(2,2) = 0.0d0

      call  Title
      if (N <= 9) then
         write(6,*) '*** Matrix ***'
         call  Show_Matrix(N, a)
      end if

      call Household_Bisec(N, a, vec, epsr, Eig, x)
      
      if (N <= 9 .and. vec) then 
        write(6,*)  '*** Mode Matrix (Eigen vectors) ***'
        call  Show_Matrix(N, x)
        do i=1, N
          write(6,*) dot_product(x(:,i),x(:,i))
        end do
      end if

      call  Show_Result(N, Eig)
      write(6,'(x,a,i5)') 'Dimension of the Matrix :', N
      write(6,*)


      stop
      end program
!====================================================================!
      subroutine  Line

      write(6,*) '----------------------------'

      return
      end subroutine
!====================================================================!
      subroutine  Title
      implicit none

      write(6,*)
      write(6,*) '**************************************************'
      write(6,*) ' Eigenvalue Problem of a Symmetric Matrix.        '
      write(6,*) ' After tridiagonalizing the matrix by Householder '
      write(6,*) ' transformation, we find engenvalues by bisection '
      write(6,*) ' method.                                          '
      write(6,*) '**************************************************'
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  B_Sort(N, t)
      implicit none
      integer, intent(in) :: N
      integer :: i, k, j
      real(8), intent(inout) :: t(N)
      real(8) :: tmax

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
      implicit none
      integer, intent(in) :: N, k
      real(8), parameter :: PI = 3.141592653589793d0
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
      character(len=30) :: FM='(2x,f12.6,a,1pd13.6)'

      call  B_Sort(N, t)
      write(6,*) '*** Results ***'
      write(6,*) '   Eigenvalue  |    Error'
      write(6,*) '-------------------------------'
      do i=1, N
          write(6,FM) t(i),'  | ', abs(SOL(N,i) - t(i))
      end do
      write(6,*)

      return
      end subroutine
