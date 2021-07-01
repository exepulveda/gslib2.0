      subroutine rand(seed,n,vector)
!---------------------------------------------------------------------
!
! This random number generator generates random numbers in ]0,1[
! Note that if the seed value is zero on the first call, a default
! value of 7931 will be used  in a linear congruential generator to
! generate 55 odd integers for the array 'itab()'. These values are
! preserved by a common statement, so that they may be used in sub-
! sequent calls by setting the seed to zero.If the value of 'seed'
! is greater than zero in a call to the subroutine, then the array
! 'itab' will be initialized and a new seed value will be returned
! by the subroutine. Best results are obtained by making the initial
! call with a seed of your choice and then setting the seed to '0'
! for all subsequent calls.
!
!---------------------------------------------------------------------
      dimension vector(*)
      common /unusual/itab(55),n1,n2,nseed
      integer rn1,seed
!
! Test to see if 55 odd integers must be generated.
!
      if((seed.gt.0).or.(nseed.lt.1)) then
            nseed = seed
            if(seed.le.0) nseed  = 7931
            do i=1,55
                  rn1=mod(nseed*9069,32768)
                  if(mod(rn1,2).eq.0) rn1 = rn1-1
                  itab(i) = rn1
                  nseed = rn1
            end do
            n1 = 0
            n2 = 24
      endif
!
! Generate "n" random components for the vector "VECTOR"
!
      do i=1,n
            itab(55-n1) = mod(itab(55-n2)*itab(55-n1),32768)
            vector(i)   = abs(float(itab(55-n1))/float(32768))
            n1 = mod(n1+1,55)
            n2 = mod(n2+1,55)
      end do
      if(seed.gt.0) seed=nseed
!
! Return with the vector of random numbers:
!
      return
      end
