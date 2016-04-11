!====================================================================!
!     Title  : trapezoid_iter.f90                                    !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-2-23-Sat                                         !
!====================================================================!
      program main
      implicit none
      integer :: j,k,kmax,N
      real*8 :: a,b,h,sumf,f
      real*8 :: T,Tp,dT,dTT,epsr
      external :: f
      parameter(kmax=100)
      parameter(a=0.0d0,b=4.0d0,epsr=1.0d-15)

      h = b - a
      Tp = (h/2.0d0)*(f(a)+f(b))
      N = 1.0d0
      k=0
      dT = Tp

      do 
         k = k + 1
         N = N*2.0d0
         h = h/2.0d0
         sumf = 0.0d0
         
         do j=1,N-1,2 
            sumf = sumf + f(a+j*h)
         end do

         T = 0.5d0*Tp + h*sumf
         dTT = dT
         dT = abs(T - Tp)

         if( dT < epsr*abs(T) ) then
              write(6,*) 'The calculation was done well.'
              exit
         end if
         if( k > kmax  ) then
              write(6,*) 'The caculation did not converge.'
              exit
         end if
         if( dT > dTT  ) then 
              write(6,*) 'Error becomes large.'
              exit
         end if
         Tp = T
      end do

      write(6,*) 'k =',k
      write(6,*) 'T =',Tp

      end program 

!====================================================================!
      real*8 function f(x)
      implicit none
      real*8 x

      f = x*x*x

      end


