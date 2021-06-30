      subroutine beyond(ivtype,nccut,ccut,ccdf,ncut,cut,cdf,zmin,zmax, &
                 ltail,ltpar,middle,mpar,utail,utpar,zval,cdfval,ierr)
!-----------------------------------------------------------------------
!
!                     Go Beyond a Discrete CDF
!                     ************************
!
! This subroutine is a general purpose subroutine to interpolate within
! and extrapolate beyond discrete points on a conditional CDF.  If the
! Z value "zval" is specified then the corresponding CDF value "cdfval"
! will be computed, if the CDF value "cdfval" is specified the
! corresponding Z value "zval" will be computed.
!
!
!
! INPUT/OUTPUT VARIABLES:
!
!   ivtype           variable type (1=continuous, 0=categorical)
!   nccut            number of cutoffs defining the conditional CDF
!   ccut()           real array of the nccut cutoffs
!   ccdf()           real array of the conditional cdf values
!   ncut             number of cutoffs defining the global CDF
!   cut()            real array of the ncut cutoffs
!   cdf()            real array of the global cdf values
!
!   zmin,zmax        minimum and maximum allowable data values
!   ltail            option to handle values in lower tail
!   ltpar            parameter required for option ltail
!   middle           option to handle values in the middle
!   mpar             parameter required for option middle
!   utail            option to handle values in upper tail
!   utpar            parameter required for option utail
!
!   zval             interesting cutoff (if -1 then it is calculated)
!   cdfval           interesting CDF (if -1 then it is calculated)
!
!
!-----------------------------------------------------------------------
      parameter(EPSLON=1.0e-20,UNEST=-1.0)
      dimension ccut(nccut),ccdf(nccut),cut(1),cdf(1)
      real      utpar,mpar,ltpar,lambda
      integer   ltail,utail,middle,cclow,cchigh
!
! Check for both "zval" and "cdfval" defined or undefined:
!
      ierr  = 1
      if(zval.gt.UNEST.and.cdfval.gt.UNEST) return
      if(zval.le.UNEST.and.cdfval.le.UNEST) return
!
! Handle the case of a categorical variable:
!
      if(ivtype.eq.0) then
            cum = 0
            do i=1,nccut
                  cum = cum + ccdf(i)
                  if(cdfval.le.cum) then
                         zval = ccut(i)
                         return
                  endif
            end do
            return
      end if
!
! Figure out what part of distribution: ipart = 0 - lower tail
!                                       ipart = 1 - middle
!                                       ipart = 2 - upper tail
      ierr  = 0
      ipart = 1
      if(zval.gt.UNEST) then
            if(zval.le.ccut(1))       ipart = 0
            if(zval.ge.ccut(nccut))   ipart = 2
      else
            if(cdfval.le.ccdf(1))     ipart = 0
            if(cdfval.ge.ccdf(nccut)) ipart = 2
      endif
!
! ARE WE IN THE LOWER TAIL?
!
      if(ipart.eq.0) then
            if(ltail.eq.1) then
!
! Straight Linear Interpolation:
!
                  powr = 1.0
                  if(zval.gt.UNEST) then
                        cdfval = powint(zmin,ccut(1),0.0,ccdf(1),zval,powr)
                  else
                        zval = powint(0.0,ccdf(1),zmin,ccut(1),cdfval,powr)
                  endif
            else if(ltail.eq.2) then
!
! Power Model interpolation to lower limit "zmin"?
!
                  if(zval.gt.UNEST) then
                        cdfval = powint(zmin,ccut(1),0.0,ccdf(1),zval,ltpar)
                  else
                        powr = 1.0 / ltpar
                        zval = powint(0.0,ccdf(1),zmin,ccut(1),cdfval,powr)
                  endif
!
! Linear interpolation between the rescaled global cdf?
!
            else if(ltail.eq.3) then
                  if(zval.gt.UNEST) then
!
! Computing the cdf value. Locate the point and the class bound:
!
                        call locate(cut,ncut,1,ncut,zval,idat)
                        call locate(cut,ncut,1,ncut,ccut(1),iupp)
!
! Straight linear interpolation if no data; otherwise, linear:
!
                        if(idat.le.0.or.idat.ge.ncut.or.iupp.le.0.or.iupp.ge.ncut) then
                               cdfval = powint(zmin,cut(1),0.0,cdf(1),zval,1.)
                         else
                               temp   = powint(cut(idat),cut(idat+1),cdf(idat),cdf(idat+1),zval,1.)
                               cdfval = temp*ccdf(1)/cdf(iupp)
                         endif
                   else
!
! Computing Z value: Are there any data out in the tail?
!
                         call locate(cut,ncut,1,ncut,ccut(1),iupp)
!
! Straight linear interpolation if no data; otherwise, local linear
! interpolation:
!
                        if(iupp.le.0.or.iupp.ge.ncut) then
                              zval = powint(0.0,cdf(1),zmin,cut(1),cdfval,1.)
                        else
                              temp = cdfval*cdf(iupp)/ccdf(1)
                              call locate(cdf,ncut,1,ncut,temp,idat)
                              if(idat.le.0.or.idat.ge.ncut) then
                                    zval = powint(0.0,cdf(1),zmin,cut(1),cdfval,1.)
                              else
                                    zval = powint(cdf(idat),cdf(idat+1),cut(idat),cut(idat+1),temp,1.)
                              end if
                        endif
                  endif
            else
!
! Error situation - unacceptable option:
!
                  ierr = 2
                  return
            endif
      endif
!
! FINISHED THE LOWER TAIL,  ARE WE IN THE MIDDLE?
!
      if(ipart.eq.1) then
!
! Establish the lower and upper limits:
!
            if(zval.gt.UNEST) then
                  call locate(ccut,nccut,1,nccut,zval,cclow)
            else
                  call locate(ccdf,nccut,1,nccut,cdfval,cclow)
            endif
            cchigh = cclow + 1
            if(middle.eq.1) then
!
! Straight Linear Interpolation:
!
                  powr = 1.0
                  if(zval.gt.UNEST) then
                        cdfval = powint(ccut(cclow),ccut(cchigh),ccdf(cclow),ccdf(cchigh),zval,powr)
                  else
                        zval = powint(ccdf(cclow),ccdf(cchigh),ccut(cclow),ccut(cchigh),cdfval,powr)
                  endif
!
! Power interpolation between class bounds?
!
            else if(middle.eq.2) then
                  if(zval.gt.UNEST) then
                        cdfval = powint(ccut(cclow),ccut(cchigh),ccdf(cclow),ccdf(cchigh),zval,mpar)
                  else
                        powr = 1.0 / mpar
                        zval = powint(ccdf(cclow),ccdf(cchigh),ccut(cclow),ccut(cchigh),cdfval,powr)
                  endif
!
! Linear interpolation between the rescaled global cdf?
!
            else if(middle.eq.3) then
                  call locate(cut,ncut,1,ncut,ccut(cclow),ilow)
                  call locate(cut,ncut,1,ncut,ccut(cchigh),iupp)
                  if(cut(ilow).lt.ccut(cclow))  ilow = ilow + 1
                  if(cut(iupp).gt.ccut(cchigh)) iupp = iupp - 1
                  if(zval.gt.UNEST) then
                        call locate(cut,ncut,1,ncut,zval,idat)
!
! Straight linear interpolation if no data; otherwise, local linear
! interpolation:
!
                        if(idat.le.0.or.idat.ge.ncut.or. &
                           ilow.le.0.or.ilow.ge.ncut.or. &
                           iupp.le.0.or.iupp.ge.ncut.or. &
                           iupp.le.ilow) then
                              cdfval=powint(ccut(cclow),ccut(cchigh), &
                                     ccdf(cclow),ccdf(cchigh),zval,1.)
                        else
                              temp = powint(cut(idat),cut(idat+1), &
                                      cdf(idat),cdf(idat+1),zval,1.)
                              cdfval=powint(cdf(ilow),cdf(iupp), &
                                     ccdf(cclow),ccdf(cchigh),temp,1.)
                        endif
                  else
!
! Straight linear interpolation if no data; otherwise, local linear
! interpolation:
!
                        if(ilow.le.0.or.ilow.ge.ncut.or. &
                           iupp.le.0.or.iupp.ge.ncut.or. &
                           iupp.le.ilow) then
                              zval=powint(ccdf(cclow),ccdf(cchigh), &
                                    ccut(cclow),ccut(cchigh),cdfval,1.)
                        else
                              temp=powint(ccdf(cclow),ccdf(cchigh), &
                                    cdf(ilow),cdf(iupp),cdfval,1.)
                              call locate(cdf,ncut,1,ncut,temp,idat)
                              if(cut(idat).lt.ccut(cclow)) idat=idat+1
                              if(idat.le.0.or.idat.ge.ncut.or. &
                                 cut(idat+1).gt.ccut(cchigh)) then
                                    zval = powint(ccdf(cclow), &
                                           ccdf(cchigh),ccut(cclow), &
                                           ccut(cchigh),cdfval,1.)
                              else
                                    zval = powint(cdf(idat),cdf(idat+1), &
                                         cut(idat),cut(idat+1),temp,1.)
                              end if
                              zval = powint(cdf(idat),cdf(idat+1), &
                                      cut(idat),cut(idat+1),temp,1.)
                        endif
                  endif
            else
!
! Error situation - unacceptable option:
!
                  ierr = 2
                  return
            endif
      endif
!
! FINISHED THE MIDDLE,  ARE WE IN THE UPPER TAIL?
!
      if(ipart.eq.2) then
            if(utail.eq.1) then
                  powr = 1.0
                  if(zval.gt.UNEST) then
                        cdfval = powint(ccut(nccut),zmax,ccdf(nccut), &
                                       1.0,zval,powr)
                  else
                        zval   = powint(ccdf(nccut),1.0,ccut(nccut), &
                                       zmax,cdfval,powr)
                  endif

            else if(utail.eq.2) then
!
! Power interpolation to upper limit "utpar"?
!
                  if(zval.gt.UNEST) then
                        cdfval = powint(ccut(nccut),zmax,ccdf(nccut), &
                                        1.0,zval,utpar)
                  else
                        powr = 1.0 / utpar
                        zval   = powint(ccdf(nccut),1.0,ccut(nccut), &
                                        zmax,cdfval,powr)
                  endif
!
! Linear interpolation between the rescaled global cdf?
!
            else if(utail.eq.3) then
                  if(zval.gt.UNEST) then
!
! Approximately Locate the point and the class bound:
!
                        call locate(cut,ncut,1,ncut,zval,idat)
                        call locate(cut,ncut,1,ncut,ccut(nccut),ilow)
                        if(cut(idat).lt.zval)        idat = idat + 1
                        if(cut(ilow).lt.ccut(nccut)) ilow = ilow + 1
!
! Straight linear interpolation if no data; otherwise, local linear
! interpolation:
!
                        if(idat.le.0.or.idat.ge.ncut.or. &
                           ilow.le.0.or.ilow.ge.ncut) then
                              cdfval = powint(ccut(nccut),zmax, &
                                        ccdf(nccut),1.0,zval,1.)
                        else
                              temp   = powint(cut(idat),cut(idat+1), &
                                        cdf(idat),cdf(idat+1),zval,1.)
                              cdfval = powint(cdf(ilow),1.0, &
                                        ccdf(nccut),1.0,temp,1.)
                        endif
                  else
!
! Computing Z value: Are there any data out in the tail?
!
                        call locate(cut,ncut,1,ncut,ccut(nccut),ilow)
                        if(cut(ilow).lt.ccut(nccut)) ilow = ilow + 1
!
! Straight linear interpolation if no data; otherwise, local linear
! interpolation:
!
                        if(ilow.le.0.or.ilow.ge.ncut) then
                              zval   = powint(ccdf(nccut),1.0, &
                                        ccut(nccut),zmax,cdfval,1.)
                        else
                              temp = powint(ccdf(nccut),1.0, &
                                      cdf(ilow),1.0,cdfval,1.)
                              call locate(cdf,ncut,1,ncut,temp,idat)
                              if(cut(idat).lt.ccut(nccut)) idat=idat+1
                              if(idat.ge.ncut) then
                                    zval   = powint(ccdf(nccut),1.0, &
                                             ccut(nccut),zmax,cdfval,1.)
                              else
                                    zval = powint(cdf(idat),cdf(idat+1), &
                                         cut(idat),cut(idat+1),temp,1.)
                              endif
                        endif
                  endif
!
! Fit a Hyperbolic Distribution?
!
            else if(utail.eq.4) then
!
! Figure out "lambda" and required info:
!
                  lambda = (ccut(nccut)**utpar)*(1.0-ccdf(nccut))
                  if(zval.gt.UNEST) then
                        cdfval = 1.0 - (lambda/(zval**utpar))
                  else
                        zval = (lambda/(1.0-cdfval))**(1.0/utpar)
                  endif
            else
!
! Error situation - unacceptable option:
!
                  ierr = 2
                  return
            endif
      endif
      if(zval.lt.zmin) zval = zmin
      if(zval.gt.zmax) zval = zmax
!
! All finished - return:
!
      return
      end
