!====================================================================!
!     Title  : Simpson_iter.f90                                      !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-2-23-Sat                                         !
!     Last modified : 2008-3-1-Sat                                   !
!                                                                    !     
!     A program which calculates an integral of a function f(x) by   !
!     using Simpson's formula. Calculation is iterated until the     !
!     value of calculation converges. After each cycle, the interval !
!     is narrowed according to h = h/2.                              !
!====================================================================!
      program main
      implicit none
      integer :: j,k,kmax,N
      real*8 :: a,b,h
      real*8 :: s1,s2,s4,ds,S,SS,epsr
      real*8,external :: f
      parameter(epsr=1.0d-15)
      parameter(kmax=100)
      parameter(a=-1.0d0,b=1.0d0)

      h = b - a
      s1 = f(a) + f(b) 
      s2 = 0.0d0
      s4 = 0.0d0

      S = (h/2.0d0)*s1
      
      k = 0
      N = 1
      
      do 
        k = k + 1
        N = N*2
        h = h/2.0d0
        s2 = s2 + s4
        s4 = 0.0d0       
      
        do j=1,N-1,2
          s4 = s4 + f(a+j*h)
        end do
        
        SS = S                                       ! last value
        S = (h/3.0d0)*(s1 + 4.0d0*s4 + 2.0d0*s2)
        ds = abs(S - SS)
         
        if( ds < epsr*abs(S) ) then
              write(6,*) 'The calculation was done well.'
              exit  
        end if
        if( k > kmax ) then 
              write(6,*) 'The calculation did not converge.'
              exit
        end if  

      write(6,*) 'h =',h
      end do

      write(6,'(1x,a,i4)')       'Number of iteration : k     =',k
      write(6,'(1x,a,1pd22.15)') 'Results             : S     =',S
      write(6,'(1x,a,1pd22.15)') 'Error               : Error =',ds

      end

!====================================================================!
      real*8 function f(x)
      implicit none
      real*8 x

      f = 2.0d0 / (x * x + 1.0d0) 
      
      end
