!====================================================================!
!     Title  : eps.f90                                               !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-2-19-Tue                                         !
!                                                                    !
!     Minimum tollerable absolute error(epsA) and machine            !
!     epsiron(epsM).                                                 !
!--------------------------------------------------------------------!
      program main 
      implicit none
      real*16 :: epsA,epsX,epsM
      integer :: i

      epsA = 1.0_16
      epsX = 0.5d0

      do while(epsX > 0.0d0)
         epsA = epsX 
         epsX = epsX/2.0d0
      end do

      epsM = 1.0d0
      epsX = 0.5d0
      
      do while(1.0d0+epsM/2.0d0 > 1.0d0)
         epsM = epsX
         epsX = epsX/2.0d0
      end do

      write(6,*) 'epsA =',epsA
      write(6,*) 'epsM =',epsM

      end
