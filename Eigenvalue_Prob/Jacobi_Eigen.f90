!====================================================================!
!     Title  : Jacobi_Eigen.f90                                      !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-5-12-Mon                                         !
!====================================================================!
      program  Main
      implicit none
      integer :: N, i, j, Step
      real(8), allocatable :: a(:,:), x(:,:)
      real(8), external :: SOL
      character(len=30) :: FM='(x,a,i3,a)'

      open(10, file='Matrix_Sym_5x5')
      read(10,*) N
      allocate(a(N,N), x(N,N))
      do i=1, N
          read(10,*) (a(i,j), j=1, N)
      end do
      close(10)

      call  Title
      call  Show_Matrix(N, a)
      call  Jacobi(N, a, x, Step)
      call  Diagonalized_Matrix(N, a)
      call  Exact_Solution
      write(6,FM) 'Iteration :', Step, ' times'
      write(6,*)
      call  Show_Result(N, a)
      call  Mode_Matrix(N, x)

      stop
      end program
!====================================================================!
      subroutine  Jacobi(N, a, x, Step)
      implicit none
      integer :: i, j, p, q
      integer, intent(in)  :: N
      integer, intent(out) :: Step
      real(8), intent(inout) :: a(N,N)
      real(8), intent(out) :: x(N,N)
      real(8) :: amax, Dmax
      real(8) :: alpha, c, s, t, tau, h
      real(8) :: apj, ajp, ajq, aqj, xpj, xjp, xqj, xjq
      real(8), parameter :: epsr = 1.0d-15
      
      Step = 0
      x = 0.0d0
      do i=1, N
          x(i,i) = 1.0d0
      end do

      do 
          p = 1
          q = 2
          amax = abs( a(p,q) )
          Dmax = abs( a(N,N) )
          do i=1, N-1
              do j=i+1, N
                  if (abs(a(i,j)) > amax) then
                      p = i
                      q = j
                      amax = abs(a(i,j))
                  end if
              end do
              Dmax = MAX(Dmax, a(i,i))
          end do
          Step = Step + 1
          if (amax < epsr * Dmax)  exit

          alpha = (a(q,q) - a(p,p)) / (2.0d0 * a(p,q))
          t     = sign(1.0d0,alpha) / &
                &   (abs(alpha) + sqrt(1.0d0 + alpha * alpha))
          c     = 1.0d0 / sqrt(1.0d0 + t * t)
          s     = t * c
          tau   = s / (1.0d0 + c)
          h     = t * a(p,q)
      
          a(p,p) = a(p,p) - h
          a(q,q) = a(q,q) + h
          a(p,q) = 0.0d0
 
          do j=1, p-1
              ajp = a(j,p)
              ajq = a(j,q)
              a(j,p) = ajp - s * (ajq + tau * ajp)
              a(j,q) = ajq + s * (ajp - tau * ajq)
          end do
          do j=p+1, q-1
              apj = a(p,j)
              ajq = a(j, q)
              a(p,j) = apj - s * (ajq + tau * apj)
              a(j,q) = ajq + s * (apj - tau * ajq)
          end do
          do j=q+1, n
              apj = a(p,j)
              aqj = a(q,j)
              a(p,j) = apj - s * (aqj + tau * apj)
              a(q,j) = aqj + s * (apj - tau * aqj)
          end do
          do j=1, n
              xjp = x(j,p)
              xjq = x(j,q)
              x(j,p) = xjp - s * (xjq + tau * xjp)
              x(j,q) = xjq + s * (xjp - tau * xjq)
          end do

      end do
      return 
      end subroutine
!====================================================================!
      subroutine  B_Sort(N, a)
      implicit none
      integer, intent(in) :: N
      integer :: i, k, j
      real(8), intent(inout) :: a(N,N)
      real(8) :: amax, temp

      amax = a(1,1)
      k = 1
      do j=1, N
          amax = a(j,j)
          k = j
          do i=j+1, N
              if (a(i,i) > amax) then
                  amax = a(i,i)
                  k = i
              end if
          end do
      
          if ( amax /= a(j,j) ) then
              a(k,k) = a(j,j)
              a(j,j) = amax
          end if
      end do

      return
      end subroutine
!====================================================================!
      function  SOL(N, k)  result(t)
      implicit none
      integer, intent(in) :: N, k
      real(8), parameter :: PI=3.141592653589793d0
      real(8) :: t

      t = 0.5d0 / (1.0d0 - &
       & cos((2.0d0*dble(k) - 1.0d0) * PI / (2.0d0*dble(n) + 1.0d0)))

      return
      end function
!====================================================================!
      subroutine  Title

      write(6,*) 
      write(6,*) '***************************************'
      write(6,*) ' Diagonalization of a symmetric matrix '
      write(6,*) ' by Jacobi method.                     '
      write(6,*) '***************************************'
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

      FM = '(x,a,' // CHAR(48+N) // 'f4.0,a)'
      write(6,*) '*** Matrix ***'
      do i=1, N
          write(6,FM) ' |', (a(i,j), j=1, N), ' |'
      end do
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  Diagonalized_Matrix(N, a)
      implicit none
      integer, intent(in) :: N
      integer :: i, j
      real(8), intent(inout) :: a(N,N)
      character(len=30) :: FM

      FM = '(x,a,' // CHAR(48+N) // 'f7.3,a)'
      write(6,*) '*** After the Transformation ***'
      do i=1, N
          write(6,FM) ' |', (a(i,j), j=1, N), ' |'
      end do
      write(6,*)

      return
      end subroutine
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
      subroutine  Show_Result(N, a)
      implicit none
      integer, intent(in) :: N
      integer :: i
      real(8), intent(in) :: a(N,N)
      real(8), external :: SOL
      character(len=30) :: FM='(2x,f9.6,a,1pd13.6)'

      call  B_Sort(N, a)
      write(6,*) '*** Results ***'
      write(6,*) ' Eigenvalue |   Error '
      write(6,*) '-------------------------------'
      do i=1, N
          write(6,FM) a(i,i),'  | ', abs(SOL(N,i) - a(i,i))
      end do
      write(6,*)

      return
      end subroutine
!====================================================================!
      subroutine  Mode_Matrix(N, x)
      implicit none
      integer, intent(in) :: N
      integer :: i, j
      real(8), intent(in) :: x(N,N)
      character(len=30) :: FM
      
      FM = '(x,a,' // CHAR(48+N) // 'f7.3,a)'
      write(6,*) '*** Mode Matrix (Eigen Vectors) ***'

      do i=1, N
          write(6,FM) ' |', (x(i,j), j=1, N), ' |'
      end do
      write(6,*)
      write(6,*) '(Note that these eigenvectors are normalized.)'
      write(6,*)

      return
      end subroutine
