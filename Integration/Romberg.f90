!====================================================================!
!     Title  : Romberg.f90                                           !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-2-27-Wed                                         !
!     Last modified  :2008-3-1-Sat                                   !
!                                                                    !
!     This program computes an integral of a function by Romberg's   !
!     method. This routine proceeds calculation with narrowing of    !
!     intervals and with extrapolation.                              !
!     If the difference of an extrapolated value and its previous    !
!     one becomes small enough, We stop extrapolation and            !
!     concentrate only on narrowing intervals.                       !
!====================================================================!
      program Romberg
!--------------------------------------------------------------------!
!     Main program.                                                  !
!     We calculate the integral by trapezoidal method at first.      !
!     Then, divide the intarval into half and calculate the integral !
!     and call subroutine to extrapolate and divide the intarval     !
!     again...                                                       !
!     If logical variable 'esc' becomes 'true',we terminate the      !
!     calculation.                                                   !
!--------------------------------------------------------------------!
      implicit none
      integer :: i,j,kmin,MAXn,N
      real*8  :: a,b,h,fsum
      real*8  :: err,S
      real*8, allocatable :: T(:)
      logical :: esc
      parameter(MAXn=20)
      parameter(a=0.0d0,b=1.0d0)
      
      allocate( T(0:MAXN) )
      h = b - a
      N = 1
      kmin = 0
      T(0) = 0.5d0*h*(f(a) + f(b))
      call pict1

 aaa: do i=1,MAXn
         write(6,*) '-----------------------------&
                    &-----------------------------'
         N = N*2
         h = 0.5d0*h
         fsum = 0.0d0

         do j=1,N-1,2 
           fsum = fsum + f(a+j*h)
         end do

         T(i) = 0.5d0*T(i-1) + h*fsum
         call extrapolation(T,S,err,esc)
         
         if(esc) then
           call pict2
           exit
         end if      
      end do  aaa

      contains
!====================================================================!
      subroutine extrapolation(T,S,err,esc)
!--------------------------------------------------------------------!
!     A subroutine that performs extrapolation. If the difference of !
!     extrapolated value and the value before extrapolation becomes  !
!     smaller than criterion, we stop the extrapolation of higer     !
!     order(by kmin = k + 1). After that we check the error only at  !
!     k=kmin.                                                        !
!     Since this subroutine is complicated, we added many comments.  !
!--------------------------------------------------------------------!
      implicit none
      real*8,intent(inout) :: T(0:MAXn)
      real*8,intent(out)   :: S,err
      logical,intent(out)  :: esc
      integer   :: j,k,m
      real*8    :: c,p,dT,dS
      real*8    :: epsr
      character :: fm1*10
      character :: fm2*30
      parameter(epsr=1.0d-15)
      parameter(fm1 = '(3f19.15)')
      parameter(fm2 = '(1x,"m :",i2," th ~",i2," th")')
     
      esc = .false.
      p=1.0d0
      
 bbb: do k=i-1,kmin,-1
        p    = p*4.0d0
        c    = 1.0d0/(p - 1.0d0)
        T(k) = T(k+1) + c*( T(k+1) - T(k) )
        dT   = abs(T(k) - T(k+1))

        if(kmin == 0) then            ! If we have not changed kmin.
           
           if(dT < epsr*abs(T(k+1))) then    
               kmin = k + 1           ! Restrict the extrapolation.
               S    = T(kmin)         ! '+1' is because 'k' increase 
               err  = dT              ! by one at each loop of 'i'.
           end if

        else if(k == kmin)  then   ! We can come to this branch only 
                                   ! after we changed kmin.
                                      ! After the change of kmin, check
                                      ! only for k=kmin.
           dS = abs(T(kmin) - S) ! Not to compare T(kmin) with extrapolated
                                 ! one, but with previous one(S).
           esc = ( dS < epsr*T(kmin) )
           if(esc) then      ! Compare 'upper' T and 'lower' T
                             ! in table.
              if( dS < err ) then
                 S   = T(kmin)
                 err = dS
              end if

              exit           ! Anyway, since the error is sufficiently 
                             ! small. Therefore exit.
           else 
           end if                 ! In this case, T(kmin) can be more
           
           kmin = kmin + 1        ! Since we have changed kmin before.
        endif
      end do   bbb

!   Write down the Romberg table on the display.
!   We use 'kmin-1' instead of 'kmin'. This is because we replaced
!   'kmin' by 'k + 1' while we compared 'T(k)' and 'T(k+1)' to 
!   change 'kmin'. We should display 'T(k)' too.
!   Note that if 'kmin = 0', we have not replaced it by 'k+1'.
      write(6,'(1x,a,i2)') 'i =',i 
      write(6,fm2) max(kmin - 1,0) ,i
      j = 0
      do
         if(max(kmin - 1,0) >= i-j-2) then                  
             write(6,fm1)  (T(i-m),m=j,i-max(kmin - 1,0))
             exit
         else
             write(6,fm1) (T(i-m),m=j,j+2)
             j = j + 3
         end if
      
      end do

      end subroutine
      
!====================================================================!
      real*8 function f(x)
      implicit none
      real*8 :: x

      f = 4.0d0/(1.0d0 + x*x)

      return
      end function

!====================================================================!
      subroutine pict1
      
      write(6,*)
      write(6,*) 'Division      ... i times.'
      write(6,*) 'Extrapolation ... m times.'
      write(6,*)
      write(6,*) '                |***************|'
      write(6,*) '                | Romberg Table |'
      write(6,*) '                *****************'
      write(6,*) '   m = 0, 3, 6        m = 1, 4           m = 2, 5'
      write(6,*) '-----------------------------&
                 &-----------------------------'
      write(6,'(1x,a,i2)') 'i =',0
      write(6,'(3f19.15)') T(0)
      
      end subroutine pict1

!====================================================================!
      subroutine pict2
      character :: Q*17
      parameter(Q = "(1x,a,1pd22.15,a)")

      write(6,*) '-----------------------------&
                 &-----------------------------'
      write(6,*) 
      write(6,*) '|**************************************************|'
      write(6,*) '|  The function we calculated.                     |'
      write(6,*) '|                                                  |'
      write(6,*) '|      1                                           |'
      write(6,*) '|     /                                            |'
      write(6,*) '|     [         4                                  |'
      write(6,*) '| S = |  dx -----------  =  pi                     |'
      write(6,*) '|     ]      1  +  x^2                             |'
      write(6,*) '|     /                  =  3.1415926535987932...  |'
      write(6,*) '|     0                                            |'
      write(6,*) '|                                                  |'
      write(6,*) '****************************************************'
      write(6,*)
      write(6,*)    '|**Results*********************************|' 
      write(6,*)    '|                                          |'
      write(6,Q)    '|    S     =', S  ,              '         |'
      write(6,Q)    '|    ERROR =', err,              '         |'
      write(6,*)    '|                                          |'
      write(6,*)    '********************************************'
  
      end subroutine
      
      end program
