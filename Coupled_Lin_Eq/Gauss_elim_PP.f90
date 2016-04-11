!====================================================================!
!     Title  : Gauss_elim_PP.f90                                     !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-3-31-Mon                                         !
!     Last modified : 2008-4-3-Thu                                   !
!====================================================================!
      program main
      implicit none
      integer :: i, j, k, n, m
      real*8, allocatable  :: a(:,:)

      open(10, file = 'Matrix_data', action='read')
      read(10,*) n, m
      allocate( a(n, n+m) )
      read(10,*) ((a(i,j), j=1, n+m), i=1, n)

      call Title
      call Draw_Matrix(n, m, a)
      call For_elim(n, m, a)
      call Back_sub(n, m, a)

      write(6,*) 'Solutions :'
      do k=n+1, n+m
          do j=1, n
              write(6,'(2x,f7.3)') a(j, k)
          end do
          write(6,*) "------------------"
      end do
      close(10)

      deallocate(a)
      end program
!====================================================================!
      subroutine  Par_Pivot(n, m, a, p, k)
      implicit none
      integer, intent(in)    :: n, m, k
      real*8,  intent(out)   :: p
      real*8,  intent(inout) :: a(n, n+m)
      real*8 :: w(n+m)
      integer :: j, L
      
      L = k
      p = a(k, k)
      
      do j=k+1, n
        if (abs(p) < abs(a(j, k))) then
          p = a(j ,k)
          L = j
        end if
      end do

!    Exchange
      if (L /= k) then
        w(:)    = a(k, :)
        a(k, :) = a(L, :)
        a(L, :) = w(:)
      end if

      end subroutine

!====================================================================!
      subroutine  For_elim(n, m, a)
      implicit none
      integer :: i, j, k, n, m
      real*8, intent(inout) :: a(n, n+m)
      real*8 :: p

      do k=1, n
          call  Par_Pivot(n, m, a, p, k)
          if ( p == 0 ) stop 'Division by 0.' 
          a(k, :) = a(k, :) / p
          do i=k+1, n
              do j=k+1, n+m
                  a(i, j) = a(i, j) - a(i, k) * a(k, j)
              end do
          end do
      end do

      return 
      end subroutine
!====================================================================!
      subroutine  Back_sub(n, m, a)
      implicit none
      integer :: i, j, k, n, m
      real*8, intent(inout) :: a(n, n+m)

      do k=n-1, 1, -1
          do i=k+1, n
              do j=n+1, n+m
                  a(k, j) = a(k, j) - a(k, i) * a(i, j)
              end do
          end do
      end do

      return 
      end subroutine
!====================================================================!
      subroutine  Draw_Matrix(n, m, a)
      implicit none
      integer,intent(in) :: n, m
      real*8, intent(in) :: a(n, n+m)
      integer :: i, j
      character(len=20) :: FM1, FM2, FM3
      parameter(FM1='(1x,a,3f8.3,a)')
      parameter(FM2='(1x,a,f8.3,a)')
      parameter(FM3='(1x,a,23x,a)')

      write(6,*)
      write(6,*) ' Coefficient Matrix'
      write(6,FM3) '--','--'
      do i=1, n
          write(6,FM1) '|',(a(i, j), j=1, n) ,' |'
      end do
      write(6,FM3) '--','--'
      write(6,*) 'Right Hand Side'
      do j=n+1, n+m
          do i=1, n
              write(6,FM2) '|',a(i, j) ,' |'
          end do
          write(6,*)   '-----------'
      end do
      write(6,*)
    
      return
      end subroutine

!====================================================================!
      subroutine  Title

      write(6,*) 
      write(6,*) '***********************************************'
      write(6,*) ' Gaussian elimination (with partial pivotting)'
      write(6,*) '***********************************************'

      return
      end subroutine
