!=====================================================================!
!     Title  : Tridiag_Bisec.f90                                      !
!     Author : Yusa Shusaku                                           !
!     Date   : 2008-5-8-Thu                                           !
!     Last modified : 2008-5-13-Tue                                   !
!=====================================================================!
      module  Global_Constant
      implicit none
      real(8), parameter :: epsa = 1.0d-100
      real(8), parameter :: epsr = 1.0d-15
      real(8), parameter :: PI = 3.141592653589793d0
      end module
!=====================================================================!
      program Main
      implicit none
      integer :: N, i, j
      real(8), allocatable :: T(:,:), Diag(:), Sub_Diag(:), Eig(:)
      real(8), external :: Exact_Sol
      character(len=30) :: FM='(2x,f9.6,a,1pd13.6)'

      open(10, file='Matrix_Tridiag_5x5')
      read(10,*) N
      allocate( T(N,N), Diag(N), Sub_Diag(N), Eig(N) )
      do i=1, N
        read(10,*) (T(i,j), j=1, N)
      end do
      close(10)

      do i=1, N-1
        Diag(i) = T(i,i)
        Sub_Diag(i) = T(i,i+1)
      end do
      Diag(N) = T(N,N)

      call Eigenvalue(N, Diag, Sub_Diag, Eig)

      call Title
      call Show_Matrix(N, T)
      write(6,*) ' *** Results ***'
      write(6,*) ' Eigenvalues  |     Error '
      call Line
      do i=1, N
        write(6,FM) Eig(i),'    |', Exact_Sol(N,i) - Eig(i) 
      end do
      write(6,*)

      stop
      end program
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
       function  Exact_Sol(N, k)  result(Lambda)
       use Global_Constant
       implicit none
       integer, intent(in) :: N
       integer :: k
       real(8) :: Lambda

      Lambda = 4.0d0 * cos(0.5d0 * dble(k) * PI / dble((N + 1))) &
             &       * cos(0.5d0 * dble(k) * PI / dble((N + 1))) 

      return
      end function
!====================================================================!
      subroutine  Line

      write(6,*) '----------------------------'

      return
      end subroutine
!====================================================================!
      subroutine  Title
      implicit none

      write(6,*)
      write(6,*) '*******************************************'
      write(6,*) ' Eigenvalue Problem of Tridiagonal Matrix. '
      write(6,*) ' Computation using Bisection Method.       '
      write(6,*) '*******************************************'
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  Show_Matrix(N, T)
      implicit none
      integer, intent(in) :: N
      integer :: i, j
      real(8), intent(in) :: T(N,N)
      character(len=30) :: FM

      FM = '(x,a,' // char(48+N) // 'f4.0,a)'
      write(6,*) '*** Matrix ***'
      do i=1, N
          write(6,FM) '  |', (T(i,j), j=1, N), ' |'
      end do
      write(6,*)

      write(6,*) '|** Exact Solution ****************|'
      write(6,*) '|                                  |'
      write(6,*) '|  N : Dimension of the matrix     | '
      write(6,*) '|                                  |'
      write(6,*) '|                 [    k * PI   ]  |'
      write(6,*) '|  t(k) = 4 cos^2 | ----------- |  |'
      write(6,*) '|                 [   2(n + 1)  ]  |'
      write(6,*) '|                                  |'
      write(6,*) '|  with k=1, 2, ... , N            |'
      write(6,*) '|                                  |'
      write(6,*) '***********************************|'
      write(6,*)
      
      return
      end subroutine
