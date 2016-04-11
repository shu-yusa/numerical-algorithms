!====================================================================!
!     Title  : DKA_method.f90                                        !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-2-20-Wed                                         !
!                                                                    !
!     A program which solves an algeblaic equation by using DKA      !
!     method.                                                        !
!     The equation is                                                !
!                                                                    !
!      f(z) = a(0)*z^n + a(1)*z^(n-1) + ... + a(n-1)*z + a(n) = 0    !
!                                                                    !
!     and 'n' stands for the order of the polynomial f(z).           !
!     The solutions z(i) are generally complex. We can see how       !
!     initial values, which are gained from subroutine 'Aberth',     !
!     approaches these solutions graphically.                        !
!                                                                    !
!     Optionally, you can makes 'n' datafiles after execution.       !
!====================================================================!
      program DKA
!--------------------------------------------------------------------!
!     Main program.                                                  !
!     You have to decide whether you make files or not.              !
!     You can calculate either default function or function that     !
!     you want to calculate.                                         !
!--------------------------------------------------------------------!
      implicit none 
      complex(8),allocatable :: a(:)
      complex(8),allocatable :: z(:)
      complex(8),external :: f
      character(len=40) :: fm='(1x,a,i2,a,f20.15," +",f20.16," i")'
      character :: P
      integer :: i,n
      logical :: opt
      
!     Make file or not.

      do 
        write(6,*) 'Do you make files ? (y/n)'
        read(5,'(a)') P
        opt = P.eq.'y'
        if( (P=='y').or.(P=='n') ) exit
      end do

!     Function to calculate.    
      
      call func_opt

!     Calculation.

      call check(n,a,z)
      call Aberth(n,a,z)
      call DK_method(n,a,opt,z)
      
!     Results.

      write(6,*)
      do i=1,n
        write(6,fm) 'z(',i,') =',z(i)
      end do

      contains

!====================================================================!
      subroutine func_opt
!--------------------------------------------------------------------!
!     Determine a function we use.                                   !
!--------------------------------------------------------------------!
 aaa: do 
        write(6,*) 'Do you calcutate default funtion ?(y/n)'
        read(5,'(a)') P
        write(6,*)
        
        if(P=='y') then
           n = 5     
           allocate(a(0:n),z(1:n))
           a = (/( 1.0d0, 0.0d0),(-1.2d0, 0.5d0),( 2.25d0,-1.25d0),   &
                 (4.05d0,5.15d0),(11.3d0,-2.6d0),(-17.4d0, -1.8d0) /)
           do i=0,n 
              write(6,fm) 'a(',i,') =',a(i)
           end do
           write(6,*)
           exit

        else if(P/='n') then
           cycle
        
        else
           write(6,*) 'Input the order of polynomial.'
           read(5,*) n
           
           if(n>20) then
               write(6,*) 'Too high order !'
               write(6,'(1x,a,/)') 'n must be such that n <= 20.'
               cycle
           end if

           allocate(a(0:n),z(1:n))

           write(6,*) 'Input the coefficients a(0) ~ a(n) in&
                     & the following form :'
           write(6,*) '(real part(dble),imaginary part(dble))'
           write(6,*)

           do i=0,n
             write(6,'(1x,a,i2,a)') 'a(',i,') ='
             read(5,*) a(i)
           end do
           
           exit
        end if
      end do aaa

      end subroutine
      
      end program

!====================================================================!
      subroutine check(n,a,z)
!--------------------------------------------------------------------!
!     Check 'a(k)' and 'n'.                                          !
!     If a(0) = 0, we rename a(1) as a(0) and repeat this procedure  !
!     until a(0) /= 0.                                               !
!     If n=1 after this procedure, the solution is only one and is   !
!     equal to -a(1)/a(0).                                           !
!--------------------------------------------------------------------!
      implicit none
      integer,intent(inout) :: n
      complex(8),intent(inout) :: a(0:n)
      complex(8),intent(out) :: z(1:n)
      integer :: i
      
      do while(a(0) == 0)
         n = n - 1         
         do i=0,n
            a(i) = a(i+1)
         end do
      end do
            
      if(n==1) then
           z(1) = - a(1)/a(0)
      else 
           return
      end if

      end subroutine

!====================================================================!
      subroutine Aberth(n,a,z)
!--------------------------------------------------------------------!
!     This subroutine gives initial values for DK method.            !
!     'ci' is the center of a circle on where we assign initial      !
!     values. The radius of the circle is given by 'R'.              !
!                                                                    !
!     Input  : n,a                                                   !
!     Output : z                                                     !
!--------------------------------------------------------------------!
      implicit none
      integer,intent(in) :: n
      complex(8),intent(in) :: a(0:n)
      complex(8),intent(out) :: z(1:n)
      integer :: i,k,L
      real(8) :: PI,R
      complex(8) :: ci,zc,b(0:n)
      parameter(PI=3.1415926535897932d0)
      
      ci = 2.0d0*PI/dble(n)*(0.0d0,1.0d0)
      zc = - a(1)/(a(0)*dble(n))
      b = a
         
      do L=1,n
         do k=1,n-L+1                     ! b(n-L) = f^(L)(zc)/l! 
            b(k) = b(k-1)*zc + b(k)       ! Horner method.
         end do
      end do

      R= 0.0d0

!    Aberth's initial value.

      do k=1,n
         R = max(R,(dble(n)*abs(b(k)/a(0)))**(1.0d0/dble(k))) 
      end do

      do i=1,n 
         z(i) = zc + R*cdexp(ci*(i-0.75))
      end do
    
      write(6,*) 'zc =',zc
      write(6,*) 'R =',R
      end subroutine

!====================================================================!
      subroutine DK_method(n,a,opt,z)
!--------------------------------------------------------------------!
!     DK iteration.                                                  !
!     Stating from Aberth's initial values, we seek for solutions by !
!     DK method.                                                     !
!     If "opt = .true.", we make file. Otherwise not.                !
!                                                                    !
!     Input  : n,a.opt,z                                             !
!     Output : z                                                     !
!--------------------------------------------------------------------!
      implicit none
      integer,intent(in) :: n
      complex(8),intent(in) :: a(0:n)
      complex(8),intent(inout) :: z(1:n)
      logical,intent(in) :: opt
      integer :: i,j
      real(8) :: epsa,epsr
      real(8) :: err,dz,d
      complex(8) :: fp,z0(1:n)
      complex(8),external :: f
      character(len=30) :: FNAME(20)
      parameter(epsa=1.0d-75,epsr=1.0d-15)

      FNAME = (/'z_1.dat','z_2.dat','z_3.dat','z_4.dat','z_5.dat',        &
                'z_6.dat','z_7.dat','z_8.dat','z_9.dat','z_10.dat',       &
                'z_11.dat','z_12.dat','z_13.dat','z_14.dat','z_15.dat',   &
                'z_16.dat','z_17.dat','z_18.dat','z_19.dat','z_20.dat'/)
      
      if(opt) then
        do i=1,n
           open(unit=10*i,file=FNAME(i))
           write(10*i,*) dble(z(i)),dimag(z(i))
        end do
      end if

 AAA: do 
         z0 = z
         err = 0.0d0
          
    BBB: do i=1,n
             fp = a(0)
             
             do j=1,n
               if(j==i) cycle 
               fp = fp*(z0(i) - z0(j))  
             end do
           
             z(i) = z0(i) - f(n,a,z0(i))/fp
            
             if(opt) write(10*i,*) dble(z(i)),dimag(z(i))

             dz = abs(z(i) - z0(i))
             d = epsa + epsr*(abs(z0(i)) + abs(z(i)))
             err = max(err,dz/d)
          end do BBB

          if (err < 1.0d0) exit       
      
      end do AAA
      
      if(opt) then
        do i=1,n
          close(10*i)
        enddo
      end if

      end subroutine

!====================================================================!
      complex(8) function f(n,a,z) result(b)
!--------------------------------------------------------------------!
!     Make a polynomial from a given set of coefficients a(k) by     !
!     Horner method.                                                 !
!--------------------------------------------------------------------!
      implicit none
      integer,intent(in) :: n
      complex(8),intent(in) :: a(0:n),z
      integer :: k

      b = a(0)
      
      do k=1,n
         b = b*z + a(k)
      end do
 
      end function 
