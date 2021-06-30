      subroutine ordrel(ivtype,ncut,ccdf,ccdfo,nviol,aviol,xviol)
!-----------------------------------------------------------------------
!
!                 Correct Order Relation Problems
!                 *******************************
!
! This subroutine identifies and corrects order relation problems in a
! conditional distribution known at a specified number of cutoffs.
!
!
!
! INPUT VARIABLES:
!
!   ivtype           variable type (0=categorical, 1=continuous)
!   ncut             number of cutoffs
!   ccdf(i)          input ccdf values
!
!
! OUTPUT VARIABLES:
!
!   ccdfo            corrected ccdf values
!   nviol()          number of order relation violations
!   aviol()          average magnitude of the order relation violations
!   xviol()          maximum magnitude of the order relation violations
!
!
!
! PROGRAMMING NOTES:
!
!   1. the arrays ccdf1 and ccdf2 are used for temporary storage of the
!      ccdf corrected sequentially upwards and downwards.  The program
!      execution will be stopped if the memory allocation of these two
!      arrays is not sufficient.
!   
!
!
!-----------------------------------------------------------------------
      parameter(MAXCUT=100)
      real      ccdf(*),ccdfo(*),aviol(*),xviol(*)
      real      ccdf1(MAXCUT),ccdf2(MAXCUT)
      integer   nviol(*)
!
! Make sure there is enough temporary storage: 
!
      if(ncut.gt.MAXCUT) then
            write(*,100) MAXCUT,ncut
 100        format('There is not enough temporary storage allocated' &
                ,/,'in subroutine ordrel: increase and recompile' &
                ,/,'      available = ',i3 &
                ,/,'      required  = ',i3)
            stop
      endif
!
! Make sure conditional cdf is within [0,1]:
!
      do i=1,ncut
            if(ccdf(i).lt.0.0) then
                  ccdf1(i) = 0.0
                  ccdf2(i) = 0.0
            else if(ccdf(i).gt.1.0) then
                  ccdf1(i) = 1.0
                  ccdf2(i) = 1.0
            else
                  ccdf1(i) = ccdf(i)
                  ccdf2(i) = ccdf(i)
            endif
      end do
!
! Correct sequentially up, then down, and then average:
!
      if(ivtype.eq.0) then
            sumcdf = 0.0
            do i=1,ncut
                  sumcdf = sumcdf + ccdf1(i)
            end do
            if(sumcdf.le.0.0) sumcdf = 1.0
            do i=1,ncut
                  ccdfo(i) = ccdf1(i) / sumcdf
            end do
      else
            do i=2,ncut
                  if(ccdf1(i).lt.ccdf1(i-1)) ccdf1(i) = ccdf1(i-1)
            end do
            do i=ncut-1,1,-1
                  if(ccdf2(i).gt.ccdf2(i+1)) ccdf2(i) = ccdf2(i+1)
            end do
            do i=1,ncut
                  ccdfo(i) = 0.5*(ccdf1(i)+ccdf2(i))
            end do
      end if
!
! Accumulate error statistics:
!
      do i=1,ncut
            if(ccdf(i).ne.ccdfo(i)) then
                  viol = abs(ccdf(i)-ccdfo(i))
                  nviol(i) = nviol(i) + 1
                  aviol(i) = aviol(i) + viol
                  xviol(i) = max(xviol(i),viol)
            endif
      end do
!
! Return with corrected CDF:
!
      return
      end
