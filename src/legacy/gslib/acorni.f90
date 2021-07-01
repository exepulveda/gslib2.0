     double precision function acorni(idum)
!-----------------------------------------------------------------------
!
! Fortran implementation of ACORN random number generator of order less
! than or equal to 12 (higher orders can be obtained by increasing the
! parameter value MAXORD).
!
!
! NOTES: 1. The variable idum is a dummy variable. The common block
!           IACO is used to transfer data into the function.
!
!        2. Before the first call to ACORN the common block IACO must
!           be initialised by the user, as follows. The values of
!           variables in the common block must not subsequently be
!           changed by the user.
!
!             KORDEI - order of generator required ( must be =< MAXORD)
!
!             MAXINT - modulus for generator, must be chosen small
!                      enough that 2*MAXINT does not overflow
!
!             ixv(1) - seed for random number generator
!                      require 0 < ixv(1) < MAXINT
!
!             (ixv(I+1),I=1,KORDEI)
!                    - KORDEI initial values for generator
!                      require 0 =< ixv(I+1) < MAXINT
!
!        3. After initialisation, each call to ACORN generates a single
!           random number between 0 and 1.
!
!        4. An example of suitable values for parameters is
!
!             KORDEI   = 10
!             MAXINT   = 2**30
!             ixv(1)   = an odd integer in the (approximate) range 
!                        (0.001 * MAXINT) to (0.999 * MAXINT)
!             ixv(I+1) = 0, I=1,KORDEI
!
!
!
! Author: R.S.Wikramaratna,                           Date: October 1990
!-----------------------------------------------------------------------
!     implicit double precision (a-h,o-z)
      parameter (KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
      common/iaco/ ixv(MAXOP1)
      do i=1,KORDEI
            ixv(i+1)=(ixv(i+1)+ixv(i))
            if(ixv(i+1).ge.MAXINT) ixv(i+1)=ixv(i+1)-MAXINT
      end do
      acorni=dble(ixv(KORDEI+1))/MAXINT
      return
      end
