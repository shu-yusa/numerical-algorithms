!====================================================================!
!     Title  : Hessenberg.f90                                        !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-5-20-Tue                                         !
!     Last modified : 2008-5-21-Wed                                  !
!                                                                    !
!     A program which transforms a matrix to Hessenberg matrix by    !
!     Householder transformation.                                    !
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
      integer :: i, j, N
      real(8), allocatable :: a(:,:), alpha(:), beta(:)
      real(8), allocatable :: T(:,:), Eig(:), P(:,:)

      open(10,file='Matrix_Sym_5x5')
      read(10,*) N
      allocate(a(N,N), alpha(N), beta(N), T(N,N), Eig(N), P(N,N))
      do i=1, N
          read(10,*) (a(i,j), j=1, N)
      end do
      close(10)
      T = a

      call  Title
      write(6,*) '*** Matrix ***'
      call  Show_Matrix(N, a)

      call  Householder(N, a, P)
      
      write(6,*) '*** Transformation Matrix P ***'
      call  Show_Matrix(N, P)
      write(6,*) '*** Hessenberg Matrix ***'
      call  Show_Matrix(N, a)

      end program
!====================================================================!
      subroutine Householder(N, a, QQ)
      implicit none
      integer :: i, k, j
      integer, intent(in) :: N
      real(8), intent(inout) :: a(N,N)
      real(8), intent(out)  :: QQ(N,N)
      real(8) :: PH(N,N), II(N,N), pq(N)
      real(8) :: t, s, w(N), aw, p(N), q(N), ps(N), qs(N), c2

      QQ = 0.0d0
      do i=1, N
          QQ(i,i) = 1.0d0
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
          q(k+1:N) = p(k+1:N) - 0.5d0 * c2 * w(k+1:N) &
                   & * Dot_product(p(k+1:N), w(k+1:N))
          ps(k:N) = c2 * Matmul(Transpose(a(k+1:N,k:N)), w(k+1:N))
          qs(1:k) = ps(1:k)
          qs(k+1:N) = ps(k+1:N) - 0.5d0 * c2 * w(k+1:N) &
                    & * Dot_product(ps(k+1:N), w(k+1:N))
          
          forall(i=1:N, j=k+1:N) a(i,j) = a(i,j) - q(i) * w(j)
          forall(i=k+1:N, j=k:N) a(i,j) = a(i,j) - w(i) * qs(j)

          pq(1:N) = c2 * Matmul(QQ(1:N,k+1:N), w(k+1:N))
          forall(i=1:N, j=k+1:N) QQ(i,j) = QQ(i,j) - pq(i) * w(j)
        end if
      end do AA
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
      write(6,*) '*************************************************'
      write(6,*) ' Transformation of a matrix to Hessenberg matrix '
      write(6,*) ' by Householder transformation.                  '
      write(6,*) '*************************************************'
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
