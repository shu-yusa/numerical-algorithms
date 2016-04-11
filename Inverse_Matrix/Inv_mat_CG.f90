!====================================================================!
!     Title  : CG_Sym.f90                                            !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-4-29-Tue                                         !
!====================================================================!
      program main
      implicit none
      integer :: i, j, N
      real(8), allocatable :: a(:,:), b(:), a_inv(:,:), c(:,:)
      real(8), allocatable :: x(:)
      character(30), parameter :: fm='(1x,a,i5,a,1pd22.15)'

      N = 200
      allocate( a(N,N), b(N), x(N), c(N,N), a_inv(N,N) )
      do i=1, N
        do j=1, i
          a(i,j) = dble(N - i + 1)
        end do
        do j=i+1, N
          a(i,j) = dble(N - j + 1)
        end do
      end do
      c = a

      if (N <= 9) then
        write(6,*) '** Original Matrix **'
        call Show_Matrix(N, a)
      end if
      call Inv_Mat(N, a, a_inv)
      if (N <= 9) then
        write(6,*) '** Inverse Matrix **'
        call Show_Matrix(N, a_inv)
        write(6,*) 'Check the product of them'
      end if
      a = matmul(c, a_inv)
      if (N <= 9) then
        call Show_Matrix(N, a)
      end if
      write(6,*) '** Check Diagonal element **'
      do i=1, N
        write(6,fm) 'i =',i,':', a(i,i)
      end do

      stop
      end program
!====================================================================!
      subroutine Inv_mat(N, a, a_inv)
      implicit none
      integer, intent(in) :: N
      integer :: i 
      real(8), intent(in) :: a(N,N)
      real(8), intent(out) :: a_inv(N,N)
      real(8), dimension(N) :: b, x
      real(8), allocatable :: c(:,:)

      allocate(c(N,N))
      c = a
      b = 0.0d0
      do i=1, N
        b(i) = 1.0d0
        call Conj_Grad(N, c, b, x)
        a_inv(:,i) = x(:)
        b(i) = 0.0d0
      end do
      deallocate(c)

      return
      end subroutine
!====================================================================!
      subroutine Conj_Grad(N, a, b, x)
      implicit none
      integer, intent(in) :: N
      integer :: i, Step
      real(8), intent(in) :: a(N,N), b(N)
      real(8), intent(out) :: x(N)
      real(8), dimension(N) :: r, p, q
      real(8) :: C1, C2, alpha, beta

      call  Guess_x(N, x)
      r = b - matmul(a, x)
      p = r
      C2 = dot_product(x, b + r)

      Step = 0
      do 
        Step = Step + 1
        C1 = C2
        q = matmul(a, p)
        alpha = dot_product(p,r) / dot_product(p,q)
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
      real(8), intent(out) :: x(N)

      do i=1 ,N
        x(i) = 10.0d0
      end do

      return
      end subroutine
!====================================================================!
      subroutine  Show_Matrix(N, a, b)
      implicit none
      integer, intent(in) :: N
      integer :: i, j
      real(8), intent(in) :: a(N,N), b(n)
      character(len=30) :: FM='(x,a,4f7.2,a,i1,a,f7.2,a)'
      
      do i=1, N
          write(6,FM) '|', (a(i,j), j=1, N), &
                     &' || x(',i,') | = |', b(i), ' |'
      end do
      write(6,*) 

      return
      end subroutine
