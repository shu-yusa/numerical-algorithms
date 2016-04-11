!====================================================================!
!     Title  : ICCG.f90                                              !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-5-2-Tue                                          !
!     imcomplete program
!====================================================================!
      program  Main
      implicit none
      integer :: N, i, j, Step
      real*8, allocatable :: a(:,:), b(:), x(:)

      open(10, file='Matrix_4x4_SO', action='read')
      read(10,*) N
      allocate( a(N,N), b(N), x(N) )
      do i=1, N
          read(10, *) (a(i,j), j=1, N), b(i)
      end do

      call  Title
      call  Show_Matrix(N, a, b)

      call  Conj_Grad(N, a, b, x, Step)
      
      call  Show_Result(N, a, b, x, Step)
      call  Show_Error(N, x)

      stop
      end program
!====================================================================!
      subroutine  IC_Decomp_1(N, A, L, D)
      implicit none
      integer, intent(in) :: N
      integer :: i, j, k
      real*8, intent(in)  :: A(N,N)
      real*8, intent(out) :: D(N), L(N,N)
      real*8 :: w

      L = a
      D = 0.0d0
      do k=1, N
          w = 0.0d0
          do j=1, k-1
              w = w + A(k,j) * A(k,j) * D(j)
          end do
          D(k) = A(k,k) - w
          L(k,k) = 1.0d0
          do j=k+1, N
              L(k,j) = 0.0d0
          end do
      end do

      return
      end subroutine
!====================================================================!
      subroutine  IC_Decomp_2(N, A, L, D)
      implicit none
      integer, intent(in) :: N
      integer :: i, j, k
      real*8, intent(in)  :: A(N,N)
      real*8, intent(out) :: L(N,N), D(N)
      real*8 :: U(N)

      L = 0.0d0
      D = 0.0d0
      do k=1, N
          do i=1, k-1
              U(i) = A(k,i) - L(i,i-1) * U(i-1) - &
                   & Dot_product(A(i,1:k-2), U(1:k-2))
          end do
           L(k,k-1)   = U(k-1) / D(k-1)
           L(k,1:k-2) = A(k,1:k-2)
           D(k) = A(k,k) - Dot_product(A(k,1:k-2), U(1:k-2)) - &
                & L(k,k-1) * U(k-1)
      end do

      return
      end subroutine
!====================================================================!
      subroutine  ICCG(N, a, b, x, Step )
      implicit none
      integer, intent(in) :: N
      integer :: i
      integer, intent(out) :: Step
      real*8, intent(in) :: A(N,N), b(N)
      real*8, intent(inout) :: x(N)
      real*8 :: r(N), p(N), q(N), L(N,N), D(N), U(N,N), C(N,N)
      real*8 :: C1, C2, alpha, beta

      call  IC_Decomp_1(N, A, L, D)
      where(D /= 0.0d0) D = sqrt(D)
      U = D * Transpose(L)
      C = Matmul(Transpose(U), U)
      call  Guess_x(N, x)
      r = b - Matmul(A, x)
      p = r
      C2 = Dot_Product(x, b + r)

      Step = 0
      do 
          Step = Step + 1
          C1 = C2
          q = Matmul(a, p)
          alpha = Dot_product(p, r) / Dot_Product(p, q)
          x = x + alpha * p
          if (Mod(Step, 3) == 0) then
              r = b - Matmul(a, x)
          else
              r = r - alpha * q
          end if
          C2 = Dot_Product(x, b + r)
          if (C1 >= C2) exit
          beta = - Dot_Product(r, q) / Dot_Product(p, q)
          p = r + beta * p
      end do
          
      return
      end subroutine
!====================================================================!
      subroutine  Guess_x(N, x)
      implicit none
      integer, intent(in) :: N
      integer :: i
      real*8, intent(out) :: x(N)

      do i=1 ,N
          x(i) = dble(i) * 3.2d0
      end do

      return
      end subroutine
!====================================================================!
      subroutine  Title

      write(6,*)
      write(6,*) '***************************************************'
      write(6,*) ' Conjugate  Gradient  Method (for Symmetic Matrix) '
      write(6,*) '***************************************************'
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  Show_Result(N, a, b, x, Step)
      implicit none
      integer, intent(in) :: N, Step
      integer :: i, j
      real*8, intent(in) :: a(N,N), b(N), x(N)
      character(len=30) :: FM ='(x,a,i3,a)'
      character(len=30) :: FM2='(x,a,i1,a,f9.5)'


      write(6,FM) 'Iteration :', Step, ' times'
      write(6,*) 
      do i=1, N
          write(6,FM2) 'x(',i,') =', x(i)
      end do
      write(6,*) 

      return
      end subroutine
!====================================================================!
      subroutine  Show_Matrix(N, a, b)
      implicit none
      integer, intent(in) :: N
      integer :: i, j
      real*8, intent(in) :: a(N,N), b(n)
      character(len=30) :: FM='(x,a,4f7.2,a,i1,a,f7.2,a)'
      
      do i=1, N
          write(6,FM) '|', (a(i,j), j=1, N), &
                     &' || x(',i,') | = |', b(i), ' |'
      end do
      write(6,*) 

      return
      end subroutine
!====================================================================!
      subroutine  Show_Error(N, x)
      implicit none
      integer, intent(in) :: N
      real*8, intent(in) :: x(N)
      character(len=30) :: FM ='(x,a,1pd12.5)'
      
      write(6,FM) 'Err(1) =',abs(x(1) - 1.0d0)
      write(6,FM) 'Err(2) =',abs(x(2) - 2.0d0)
      write(6,FM) 'Err(3) =',abs(x(3) + 1.0d0)
      write(6,FM) 'Err(4) =',abs(x(4) - 1.0d0)
      write(6,*)

      return
      end subroutine
