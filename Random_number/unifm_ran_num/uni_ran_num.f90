!======================================================================!
!     Title  : uni_ran_num.f90                                         !
!     Code   : Taken from "Numerical Recipes in Fortran 90"            !
!     Date   : 2009-2-9-Mon                                            !
!                                                                      !
!     This code is taken from "Numerical Recipes in Fortran 90".       !
!======================================================================!
      subroutine ran(idum, r)
!----------------------------------------------------------------------!
!     "Minimal random number generator of Park and Miller combined     !
!     with a Marsaglia shift sequence. Returns a uniform random        !
!     deviate between 0.0 and 1.0 (exclusive of the endpoint values).  !
!     This fully portable, scalar generator has the "traditional"      !
!     (not Fortran 90) calling squence with a random deviate as the    !
!     returned function value: call with idum a negative integer to    !
!     initialize; thereafter, do not alter idum except to reinitialize.!
!     The period of this generator is about 3.1 X 10^18.               !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(inout) :: idum
      integer, parameter :: IA=16807,IM=2147483647,IQ=12773,IR=2836
      integer, save :: ix=-1, iy=-1, k
      real*8, intent(out) :: r
      real*8, save :: am

      if (idum <= 0 .or. iy < 0) then
         am = nearest(1.0d0, -1.0d0) / IM
         iy = ior(ieor(888889999, abs(idum)), 1)
         ix = ieor(777755555, abs(idum))
         idum = abs(idum) + 1
      end if
      ix = ieor(ix, ishft(ix,13))
      ix = ieor(ix, ishft(ix,-17))
      ix = ieor(ix, ishft(ix,5))
      k = iy / IQ
      iy = IA * (iy - k * IQ) - IR * k
      if (iy < 0) iy = iy + IM
      r = am * ior(iand(IM, ieor(ix,iy)), 1)

      end subroutine
