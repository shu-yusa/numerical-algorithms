!====================================================================!
!     Title  : Horner_method.f90                                     !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-2-20-Wed                                         !
!                                                                    !
!     Derivative of a n-th order polynomial f(x).                    !
!     The calculation is based on Horner method.                     !
!--------------------------------------------------------------------!
      program main
      implicit none
      integer :: n,L,k
      real*8 x0
      real*8,allocatable :: a(:),b(:)
      character(len=30) :: fm
      parameter(fm='(1x,a,i2,a,f8.4)')
      parameter(n=4)
      parameter(x0=0.0d0)

      allocate(a(0:n),b(0:n))
      a = (/1.0d0,-3.0d0,5.0d0,1.0d0,2.0d0/)

      call pict
      call Horner(b)
       
      do k=0,n
        write(6,fm) 'f^(',k,')(0) = ',b(n-k)
      end do

      contains 
!====================================================================!
      subroutine Horner(b)
!--------------------------------------------------------------------!
!     Calculation of derivatives by using Horner method.             !
!     b(n-k) : N-th derivative of f(x) evaluated at x=x0             !
!--------------------------------------------------------------------!
      implicit none
      real*8,intent(out) :: b(0:n)
      integer :: k,L,fact
      
      do k=0,n
         b(k) = a(k)
      end do

      do L=1,n
         do k=1,n-L+1
            b(k) = b(k) + b(k-1)*x0
         end do
      end do

      fact = 1
      
      do k=2,n
         fact = fact*k
         b(n-k) = b(n-k)*dble(fact)   ! b(n-k) = f^(k)(x0) 
      end do

      end subroutine 

!====================================================================!
      subroutine pict

      write(6,*)
      write(6,*) 'f(x) = x^4 - 3x^3 + 5x^2 + x + 2'
      write(6,*)
      write(6,*) 'Values of derivative at x=0.'

      end subroutine

      end program
