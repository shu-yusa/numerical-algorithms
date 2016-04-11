!====================================================================!
!     Title  : CG_GL.f90                                             !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-5-2-Fri                                          !
!====================================================================!
      program  Main
      implicit none
      integer :: N, i, j, Step
      real(8), allocatable :: a(:,:), b(:), x(:)

      open(10, file='Matrix_3x3', action='read')
      read(10,*) N
      allocate( a(N,N), b(N), x(N) )
      do i=1, N
          read(10, *) (a(i,j), j=1, N), b(i)
      end do

      call  Title
      call  Show_Matrix(N, a, b)

      call  Conj_Grad(N, a, b, x, Step)
      
      call  Show_Result(N, a, b, x, Step)

      stop
      end program
!====================================================================!
      subroutine  Conj_Grad(N, a, b, x, Step )
      implicit none
      integer, intent(in) :: N
      integer :: i
      integer, intent(out) :: Step
      real(8), intent(in) :: a(N,N), b(N)
      real(8), intent(inout) :: x(N)
      real(8) :: r(N), p(N), q(N), c(N), s(N)
      real(8) :: C1, C2, alpha, beta

      call  Guess_x(N, x)
      c = b - Matmul(a, x)
      r = Matmul(Transpose(a), c)
      p = r
      C2 = Dot_Product(c, c)

      Step = 0
      do 
          Step = Step + 1
          C1 = C2
          q = Matmul(a, p)
          s = Matmul(Transpose(a), q)
          alpha = Dot_product(p, r) / Dot_Product(q, q)
          x = x + alpha * p
          c = c - alpha * q
          r = r - alpha * s
          C2 = Dot_Product(c, c)
          if (C1 <= C2) exit
          beta = - Dot_Product(r, s) / Dot_Product(q, q)
          p = r + beta * p
      end do
          
      return
      end subroutine
!====================================================================!
      subroutine  Guess_x(N, x)
      implicit none
      integer, intent(in) :: N
      integer :: i
      real(8), intent(out) :: x(N)

      do i=1 ,N
          x(i) = dble(i) * 3.2d0
      end do

      return
      end subroutine
!====================================================================!
      subroutine  Title

      write(6,*)
      write(6,*) '*******************************'
      write(6,*) '  Conjugate  Gradient  Method  '
      write(6,*) '*******************************'
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  Show_Result(N, a, b, x, Step)
      implicit none
      integer, intent(in) :: N, Step
      integer :: i, j
      real(8), intent(in) :: a(N,N), b(N), x(N)
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
      real(8), intent(in) :: a(N,N), b(n)
      character(len=30) :: FM='(x,a,3f7.2,a,i1,a,f7.2,a)'
      
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
      real(8), intent(in) :: x(N)
      character(len=30) :: FM ='(x,a,1pd12.5)'
      
      write(6,FM) 'Err(1) =',abs(x(1) - 1.0d0)
      write(6,FM) 'Err(2) =',abs(x(2) - 2.0d0)
      write(6,FM) 'Err(3) =',abs(x(3) + 1.0d0)
      write(6,FM) 'Err(4) =',abs(x(4) - 1.0d0)
      write(6,*)

      return
      end subroutine
