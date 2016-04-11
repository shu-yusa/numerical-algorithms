!====================================================================!
!     Title  : Gauss_elim_CP.f90                                     !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-3-31-Mon                                         !
!====================================================================!
      program main
      implicit none
      integer :: i, j, k, n, m
      integer, allocatable :: mm(:)
      real(8), allocatable  :: a(:,:)

      open(10, file = 'Matrix_data', action='read')
      read(10,*) n, m
      allocate( a(n, n+m), mm(1:n) )
      read(10,*) ((a(i,j), j=1, n+m), i=1, n)

      call Title
      call Draw_Matrix(n, m, a)

      call For_elim(n, m, a, mm)
      call Back_sub(n, m, a, mm)

      write(6,*) 'Solutions :'
      do k=n+1, n+m
          do j=1, n
              write(6,'(2x,f7.3)') a(j, k)
          end do
      end do

      close(10)
      deallocate(a)
      deallocate(mm)
      end program
!====================================================================!
      subroutine  Cmp_Pivot(n, m, a, mm, p, k)
      implicit none
      integer, intent(in)    :: n, m, k
      integer, intent(inout) :: mm(1:n)
      integer :: i, j, IR, IC, iw
      real(8), intent(out)   :: p
      real(8), intent(inout) :: a(n,n+m)
      real(8)  :: w(n+m)

      IR = k
      IC = k
      p = a(k,k)
 
      do i=k, n
        do j=k, n
          if (abs(p) < abs(a(i,j))) then
            IR = i
            IC = j
            p = a(i,j)
          end if
        end do
      end do

!    Exchange
      if (IR /= k) then   ! Row
        w(k:n+m) = a(IR,k:n+m)
        a(IR,k:n+m) = a(k,k:n+m)
        a(k,k:n+m) = w(k:n+m)
      end if
      if (IC /= k) then   ! Column
        iw = mm(k)
        mm(k) = mm(IC)
        mm(IC) = iw
        w(1:n) = a(1:n,IC)
        a(1:n,IC) = a(1:n,k)
        a(1:n,k) = w(1:n)
      end if

      end subroutine
!====================================================================!
      subroutine  swap_real(x, y)
      implicit none
      real(8), intent(inout) :: x, y
      real(8) :: w

      w = x
      x = y
      y = w

      return
      end subroutine 
!====================================================================!
      subroutine  swap_int(p, q)
      implicit none
      integer, intent(inout) :: p, q
      integer :: r

      r = p
      p = q
      q = r

      return
      end subroutine
!====================================================================!
      subroutine  For_elim(n, m, a, mm)
      implicit none
      integer :: i, j, k, n, m
      integer, intent(out)  :: mm(1:n)
      real(8), intent(inout) :: a(n, n+m)
      real(8) :: p

      do k=1, n
        mm(k) = k
      end do

      do k=1, n
        call  Cmp_Pivot(n, m, a, mm, p, k)
        if ( p == 0 ) stop 'Division by 0.' 
        a(k, k+1:n+m) = a(k, k+1:n+m) / p
        do i=k+1, n
          a(i,k+1:n+m) = a(i,k+1:n+m) - a(i,k) * a(k,k+1:n+m)
        end do
      end do

      return 
      end subroutine
!====================================================================!
      subroutine  Back_sub(n, m, a, mm)
      implicit none
      integer, intent(in) :: n, m
      integer, intent(inout) :: mm(n)
      integer :: i, j, k, ik, iw
      real(8), intent(inout) :: a(n, n+m)
      real(8) :: w(n+m)

      do k=n-1, 1, -1
        a(k,n+1:n+m) = a(k,n+1:n+m)- matmul(a(k,k+1:n),a(k+1:n,n+1:n+m))
      end do

      do k=1, n-1
        do while(mm(k) /= k)
          ik = mm(k)
          mm(k) = mm(ik)
          mm(ik) = ik
          w = a(ik,:)
          a(ik,:) = a(k,:)
          a(k,:) = w
        end do
      end do

      return 
      end subroutine
!====================================================================!
      subroutine  Draw_Matrix(n, m, a)
      implicit none
      integer,intent(in) :: n, m
      real(8), intent(in) :: a(n, n+m)
      integer :: i, j
      character :: FM1*20, FM2*20
      parameter(FM1='(1x,a,3f8.3,a)')
      parameter(FM2='(1x,a,f8.3,a)')

      write(6,*)
      write(6,*) ' Coefficient Matrix'
!     write(6,'(1x,a,23x,a)') '--','--'
      do i=1, n
          write(6,FM1) '|',(a(i, j), j=1, n+m-1) ,' |'
      end do
      write(6,'(1x,a,23x,a)') '--','--'
      write(6,*) 'Right Hand Side'
      do i=1, n
          write(6,FM2) '|',a(i, n+1) ,' |'
      end do
      write(6,*)

    

      return
      end subroutine

!====================================================================!
      subroutine  Title
      implicit none

      write(6,*)
      write(6,*) '************************************************'
      write(6,*) ' Gaussian elimination (with complete pivotting)'
      write(6,*) '************************************************'
      
      return 
      end subroutine 
