      subroutine nscore(nd,vr,tmin,tmax,iwt,wt,tmp,lout,vrg,ierror)
!-----------------------------------------------------------------------
!
!              Transform Univariate Data to Normal Scores
!              ******************************************
!
! This subroutibe takes "nd" data "vr(i),i=1,...,nd" possibly weighted
! by "wt(i),i=,...,nd" and returns the normal scores transform N(0,1)
! as "vrg(i),i=1,...,nd".  The extra storage array "tmp" is required
! so that the data can be returned in the same order (just in case
! there are associated arrays like the coordinate location).
!
!
!
! INPUT VARIABLES:
!
!   nd               Number of data (no missing values)
!   vr(nd)           Data values to be transformed
!   tmin,tmax        data trimming limits
!   iwt              =0, equal weighted; =1, then apply weight
!   wt(nd)           Weight for each data (don't have to sum to 1.0)
!   tmp(nd)          Temporary storage space for sorting
!   lout             if > 0 then transformation table will be written
!
!
!
! OUTPUT VARIABLES:
!
!   vrg(nd)          normal scores
!   ierror           error flag (0=error free,1=problem)
!
!
!
! EXTERNAL REFERENCES:
!
!   gauinv           Calculates the inverse of a Gaussian cdf
!   sortem           sorts a number of arrays according to a key array
!
!
!
!-----------------------------------------------------------------------
      parameter(EPSLON=1.0e-20)
      real      vr(nd),wt(nd),vrg(nd),tmp(nd)
      real*8    pd
!
! Sort the data in ascending order and calculate total weight:
!
      ierror = 0
      twt    = 0.0
      do i=1,nd
            tmp(i) = real(i)
            if(vr(i).ge.tmin.and.vr(i).lt.tmax) then
                  if(iwt.eq.0) then
                        twt = twt + 1.
                  else
                        twt = twt + wt(i)
                  end if
            end if
      end do
      if(nd.lt.1.or.twt.lt.EPSLON) then
            ierror = 1
            return
      end if
      call sortem(1,nd,vr,2,wt,tmp,d,e,f,g,h)
!
! Compute the cumulative probabilities:
!
      oldcp = 0.0
      cp    = 0.0
      do i=1,nd
            cp     =  cp + wt(i) / twt
            wt(i)  = (cp + oldcp)/ 2.0
            oldcp  =  cp
            call gauinv(dble(wt(i)),vrg(i),ierr)
            if(lout.gt.0) write(lout,'(f12.5,1x,f12.5)') vr(i),vrg(i)
      end do
!
! Get the arrays back in original order:
!
      call sortem(1,nd,tmp,3,wt,vr,vrg,e,f,g,h)
!
! Finished:
!
      return
      end
