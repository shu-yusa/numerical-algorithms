!====================================================================!
!     Title  : Newton_method.f90                                     !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-2-19-Tue                                         !
!                                                                    !
!     Seek for a solution of an following eqation by Newton method : !
!         (x-1)^2*x=0                                                !
!====================================================================!
      program main
      implicit none
      
      write(6,*) 'Convergence of order one'
      call Newton
      write(6,*)  
      write(6,*) 'Convergence of order two'
      call Newton_im

      end 

!====================================================================!
      subroutine Newton
!--------------------------------------------------------------------!
!     Ordinary Newton method.                                        !
!--------------------------------------------------------------------!

      implicit none
      integer :: i
      real*8 :: epsa,epsr,err
      real*8 :: x,x0,df
      real*8,external :: f
      character :: fm*30 
      parameter(epsa=1.0d-75,epsr=1.0d-11)
      parameter(fm = '(1x,a,i2,a,f13.10)')

      x= 1.3d0
      i=0
      write(6,fm) 'x(',i,') =',x

      do
         i = i + 1
         x0 = x
         x = x0 - f(x0,df)/df
         err = abs(x - x0) / (epsa + epsr * (abs(x) + abs(x0)))
         write(6,fm) 'x(',i,') =',x
         if(err < 1.0d0) exit
      
      end do

      end subroutine

!====================================================================!
      subroutine Newton_im
!--------------------------------------------------------------------!
!     Newton_method for an equation which has 'juukai' of order      !
!     two. (I do not know the corresponding word in English)         !
!--------------------------------------------------------------------!
      implicit none
      integer :: i
      real*8 :: epsa,epsr,err
      real*8 :: x,x0,df
      real*8,external :: f
      character :: fm*30 
      parameter(epsa=1.0d-75,epsr=1.0d-11)
      parameter(fm = '(1x,a,i2,a,f13.10)')

      x= 1.3d0
      i=0
      write(6,fm) 'x(',i,') =',x

      do
         i = i + 1
         x0 = x
         x = x0 - 2.0d0*f(x0,df)/df
         err = abs(x - x0) / (epsa + epsr * (abs(x) + abs(x0)))
         write(6,fm) 'x(',i,') =',x
         if(err < 1.0d0) exit
      
      end do

      end subroutine

!====================================================================!
      real*8 function f(x,df)
!--------------------------------------------------------------------!
!      Definitions of a function and its derivative.                 !
!--------------------------------------------------------------------!
      implicit none
      real*8,intent(in) :: x
      real*8,intent(out) :: df

      f = (x - 1.0d0) ** 2 * x
      df = 3.0d0 * ((x - 4.0d0 / 3.0d0) * x + 1.0d0 / 3.0d0)

      end function

