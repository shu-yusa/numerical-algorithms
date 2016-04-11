!====================================================================!
!     Title  : Adms_Bf_varh.f90                                      !
!     Author : Yusa Shusaku                                          !
!     Date   : 2008-3-30-Sun                                         !
!====================================================================!
      module  Com_var
      implicit none
      end module
!====================================================================!
      program main
      implicit none

      end program
!====================================================================!
      subroutine  Machine_eps
      implicit none
      real*8 :: epsm = 0.03125d0
      real*8 :: epsm4

      do while(1.0d0 + epsm > 1.0d0)
          epsm  = epsm * 2.0d0
          epsm4 = epsm * 4.0d0
      end do
      
      end subroutine
!====================================================================!
      subroutine  
