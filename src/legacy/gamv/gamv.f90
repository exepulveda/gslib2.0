!
! Module to declare dynamic arrays in multiple subroutines:
!
      module geostat_gamv
      
      real,allocatable         :: x(:),y(:),z(:),vr(:,:),azm(:),atol(:), &
                                  bandwh(:),dip(:),dtol(:),bandwd(:)
      real*8,allocatable       :: sills(:),dis(:),gam(:),hm(:), &
                                  tm(:),hv(:),tv(:),np(:)
      integer,allocatable      :: ivtail(:),ivhead(:),ivtype(:)
      character*12,allocatable :: names(:)

      real      EPSLON,VERSION
      real      xlag,xltol,tmin,tmax
      integer   nd,nlag,ndir,nvarg,isill,test
      character outfl*512
      
      end module
!
!
!
      program main
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %
! Junior University.  All rights reserved.                             %
!                                                                      %
! The programs in GSLIB are distributed in the hope that they will be  %
! useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
! responsibility to anyone for the consequences of using them or for   %
! whether they serve any particular purpose or work at all, unless he  %
! says so in writing.  Everyone is granted permission to copy, modify  %
! and redistribute the programs in GSLIB, but only under the condition %
! that this notice and the above copyright notice remain intact.       %
!                                                                      %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
!               Variogram of Irregularly Spaced 3-D Data
!               ****************************************
!
! This is a template driver program for GSLIB's "gamv" subroutine. The
! input data must be entered with coordinates in a GEOEAS format file.
! The User's Guide can be referred to for more details.
!
! The program is executed with no command line arguments.  The user
! will be prompted for the name of a parameter file.  The parameter
! file is described in the documentation (see the example gamv.par)
!
!
!
! The output file will contain each directional variogram ordered by
! direction and then variogram (the directions cycle fastest then the
! variogram number).  For each variogram there will be a one line
! description and then "nlag" lines with the following:
!
!        a) lag number (increasing from 1 to nlag)
!        b) separation distance
!        c) the "variogram" value
!        d) the number of pairs for the lag
!        e) the mean of the data contributing to the tail
!        f) the mean of the data contributing to the head
!
!
!
! AUTHOR: Clayton V. Deutsch                             DATE: 1989-1999
!-----------------------------------------------------------------------
      use geostat_gamv
      EPSLON  = 1.0e-20
      VERSION = 2.905
!
! Read the Parameter File:
!
      call readparm
!
! Call gamv to compute the required variograms:
!
      call gamv
!
! Write Results:
!
      call writeout
!
! Finished:
!
      write(*,9998) VERSION
 9998 format(/' GAMV Version: ',f5.3, ' Finished'/)
      stop
      end
 
 
 
      subroutine readparm
!-----------------------------------------------------------------------
!
!                  Initialization and Read Parameters
!                  **********************************
!
! The input parameters and data are read in from their files. Some quick
! error checking is performed and the statistics of all the variables
! being considered are written to standard output.
!
!
!
!-----------------------------------------------------------------------
      ! use msflib
      use geostat_gamv
      parameter(MV=500)

      real      var(MV),cut(MV)
      real*8    avg(MV),ssq(MV)
      integer   ivar(MV),num(MV),ivc(MV),indflag(MV)
      character datafl*512,str*512
      logical   testfl,testdat
      real,allocatable :: vrmin(:),vrmax(:)

      data      lin/1/,ncut/0/
!
! Note VERSION number:
!
      write(*,9999) VERSION
 9999 format(/' GAMV Version: ',f5.3/)
!
! Get the name of the parameter file - try the default name if no input:
!
      do i=1,512
            str(i:i) = ' '
      end do
      call getarg(1,str)
      if(str(1:1).eq.' ')then
            write(*,*) 'Which parameter file do you want to use?'
            read (*,'(a)') str
      end if
      if(str(1:1).eq.' ') str(1:20) = 'gamv.par            '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'gamv.par            ') then
                  write(*,*) '        creating a blank parameter file'
                  call makepar
                  write(*,*)
            end if
            stop
      endif
      open(lin,file=str,status='OLD')
!
! Find Start of Parameters:
!
 1    read(lin,'(a4)',end=98) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
!
! Read Input Parameters:
!
      read(lin,'(a512)',err=98) datafl
      call chknam(datafl,512)
      write(*,*) ' data file = ',datafl(1:40)

      read(lin,*,err=98) ixl,iyl,izl
      write(*,*) ' columns for X,Y,Z = ',ixl,iyl,izl

      read(lin,*,err=98) nvar
      write(*,*) ' number of variables = ',nvar
      backspace lin

      read(lin,*,err=98) j,(ivar(i),i=1,nvar)
      write(*,*) ' columns = ',(ivar(i),i=1,nvar)

      read(lin,*,err=98) tmin,tmax
      write(*,*) ' trimming limits = ',tmin,tmax

      read(lin,'(a512)',err=98) outfl
      call chknam(outfl,512)
      write(*,*) ' output file = ',outfl(1:40)

      read(lin,*,err=98) nlag
      write(*,*) ' number of lags = ',nlag
      if(nlag.lt.1)       stop 'nlag is too small: check parameters'

      read(lin,*,err=98) xlag
      write(*,*) ' lag distance = ',xlag
      if(xlag.le.0.0) stop 'xlag is too small: check parameter file'

      read(lin,*,err=98) xltol
      write(*,*) ' lag tolerance = ',xltol

      read(lin,*,err=98) ndir
      write(*,*) ' number of directions = ',ndir
      
      allocate (azm(ndir),stat = test)   
      allocate (atol(ndir),stat = test)
      allocate (bandwh(ndir),stat = test)
      allocate (dip(ndir),stat = test)
      allocate (dtol(ndir),stat = test)
      allocate (bandwd(ndir),stat = test)      

      if(ndir.lt.1)       stop 'ndir is too small: check parameters'

      do i=1,ndir
            read(lin,*,err=98) azm(i),atol(i),bandwh(i), &
                               dip(i),dtol(i),bandwd(i)
            write(*,*) ' azm, atol, bandwh = ',azm(i),atol(i),bandwh(i)
            write(*,*) ' dip, dtol, bandwd = ',dip(i),dtol(i),bandwd(i)
            if(bandwh(i).lt.0.0) then
                  write(*,*) ' Horizontal bandwidth is too small!'
                  stop
            endif
            if(bandwd(i).lt.0.0) then
                  write(*,*) ' Vertical bandwidth is too small!'
                  stop
            endif
      end do

      read(lin,*,err=98) isill
      write(*,*) ' flag to standardize sills = ',isill

      read(lin,*,err=98) nvarg
      write(*,*) ' number of variograms = ',nvarg
      if(nvarg.lt.1)      stop 'nvarg is too small: check parameters'      

      mxdlv=ndir*(nlag+2)*nvarg
      allocate (dis(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if      
      allocate (gam(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if      
      allocate (hm(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if      
      allocate (tm(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if      
      allocate (hv(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if      
      allocate (tv(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if      
      allocate (np(mxdlv),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if      
      allocate (ivtail(nvarg),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if      
      allocate (ivhead(nvarg),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if      
      allocate (ivtype(nvarg),stat=test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if      
      allocate (names(nvar+nvarg),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if
      
      ncut = 0
      do i=1,nvarg
            read(lin,*,err=98) ivtail(i),ivhead(i),ivtype(i)
            write(*,*) ' tail,head,type = ', &
                         ivtail(i),ivhead(i),ivtype(i)
            if(ivtype(i).eq.9.or.ivtype(i).eq.10) then
                   ncut = ncut + 1
                   if(tmin.gt.0.0) stop 'tmin interferes with indicators!'
                   if(tmax.le.1.0) stop 'tmax interferes with indicators!'
                   backspace lin
                   read(lin,*,err=98) ii,jj,kk,cut(ncut)
                   if(ivtype(i).eq.9)  indflag(ncut) = 1
                   if(ivtype(i).eq.10) indflag(ncut) = 0
                   ivc(ncut) = ivtail(i)
                   ivtail(i) = nvar + ncut
                   ivhead(i) = nvar + ncut
                   write(names(nvar+ncut),140) ncut
 140               format('Indicator ',i2)
                   write(*,*) ' indicator threshold = ',cut(ncut)
            endif
      end do
      write(*,*)
      close(lin)
      MAXVAR = nvar + ncut
!
! Perform some quick error checking:
!
      if(xltol.le.0.0) then
            write(*,*) 'xltol is too small: resetting to xlag/2'
            xltol = 0.5*xlag
      endif
!
! Check to make sure the data file exists, then either read in the
! data or write an error message and stop:
!
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR data file ',datafl,' does not exist!'
            stop
      endif
!
! The data file exists so open the file and read in the header
! information. Initialize the storage that will be used to summarize
! the data found in the file:
!
      open(lin,file=datafl,status='OLD')

      read(lin,*,err=99)
      read(lin,*,err=99)       nvari
      do i=1,nvari
            read(lin,*)
      end do
      maxdat = 0
 22   read(lin,*,end=44,err=99)(var(j),j=1,nvari)
      maxdat = maxdat + 1
      go to 22
 44   continue
      write(*,*)'maxdat = ',maxdat

      allocate (vr(maxdat,MAXVAR),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if

      allocate (vrmin(MAXVAR),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if

      allocate (vrmax(MAXVAR),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if

      allocate (sills(MAXVAR),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if

      allocate (x(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if

      allocate (y(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if

      allocate (z(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if

      rewind(lin) 
      read(lin,'(a40)',err=99) str
      read(lin,*,err=99)       nvari
      do i=1,nvari
            read(lin,'(a40)',err=99) str
            do iv=1,nvar
                  j=ivar(iv)
                  if(i.eq.j) names(iv) = str(1:12)
            end do
            num(i) = 0
            avg(i) = 0.0
            ssq(i) = 0.0
      end do
!
! Read all the data until the end of the file:
!
      nd = 0
 2    continue
      read(lin,*,end=9,err=99) (var(j),j=1,nvari)
      testdat = .false.
      do iv=1,nvar
            j=ivar(iv)
            if(var(j).ge.tmin.and.var(j).lt.tmax) testdat = .true.
      end do
      if(.not.testdat) go to 2
      nd = nd + 1
!
! Acceptable data, make sure there are not too many data:
!
      do iv=1,nvar
            j=ivar(iv)
            vr(nd,iv) = var(j)
            if(var(j).ge.tmin.and.var(j).lt.tmax) then
                  num(iv) = num(iv) + 1
                  avg(iv) = avg(iv) + dble(var(j))
                  ssq(iv) = ssq(iv) + dble(var(j)*var(j))
            endif
      end do
      if(ixl.le.0) then
            x(nd) = 0.0
      else
            x(nd) = var(ixl)
      endif
      if(iyl.le.0) then
            y(nd) = 0.0
      else
            y(nd) = var(iyl)
      endif
      if(izl.le.0) then
            z(nd) = 0.0
      else
            z(nd) = var(izl)
      endif
      go to 2
 9    continue
      close(lin)
!
! Compute the averages and variances as an error check for the user:
!
      do iv=1,nvar
            sills(iv) = -999.
            if(num(iv).gt.0) then
                  avg(iv)   = avg(iv)/dble(num(iv))
                  ssq(iv)   =(ssq(iv)/dble(num(iv)))-avg(iv)*avg(iv)
                  sills(iv) = ssq(iv)
                  write(*,*) 'Variable number ',iv
                  write(*,*) '  Number   = ',num(iv)
                  write(*,*) '  Average  = ',real(avg(iv))
                  write(*,*) '  Variance = ',real(ssq(iv))
            endif
      end do
!
! Construct Indicator Variables if necessary:
!
      do ic=1,ncut
            iv   = ivc(ic)
            jv   = nvar + ic
            ptot = 0.0
            p1   = 0.0
            do id=1,nd
                  if(vr(id,iv).le.tmin.or.vr(id,iv).gt.tmax) then
                        vr(id,jv) = tmin - EPSLON
                  else
                        if(indflag(ic).eq.1) then
                              if(vr(id,iv).lt.cut(ic)) then
                                    vr(id,jv) = 0.0
                              else
                                    vr(id,jv) = 1.0
                              endif
                              p1   = p1   + vr(id,jv)
                              ptot = ptot + 1.0
                        else
                              vr(id,jv) = 0.0
                              if(int(vr(id,iv)+0.5).eq.int(cut(ic)+0.5)) &
                              vr(id,jv) = 1.0
                              p1   = p1   + vr(id,jv)
                              ptot = ptot + 1.0
                        end if
                  end if
            end do
            p1        = p1 / max(ptot,1.0)
            sills(jv) = dble (p1*(1.0-p1))
      end do
!
! Establish minimums and maximums:
!
      do i=1,MAXVAR
            vrmin(i) =  1.0e21
            vrmax(i) = -1.0e21
      end do
      do id=1,nd
            do iv=1,nvar+ncut
                  if(vr(id,iv).ge.tmin.and.vr(id,iv).lt.tmax) then
                        if(vr(id,iv).lt.vrmin(iv)) vrmin(iv) = vr(id,iv)
                        if(vr(id,iv).gt.vrmax(iv)) vrmax(iv) = vr(id,iv)
                  end if
            end do
      end do
!
! Check on the variogams that were requested:
!
      call check(vrmin,vrmax)
      return
!
! Error in an Input File Somewhere:
!
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'
      end



      subroutine gamv
!-----------------------------------------------------------------------
!
!              Variogram of 3-D Irregularly Spaced Data
!              ****************************************
!
! This subroutine computes a variety of spatial continuity measures of a
! set for irregularly spaced data.  The user can specify any combination
! of direct and cross variograms using any of eight "variogram" measures
!
!
!
! INPUT VARIABLES:
!
!   nd               Number of data (no missing values)
!   x(nd)            X coordinates of the data
!   y(nd)            Y coordinates of the data
!   z(nd)            Z coordinates of the data
!   nv               The number of variables
!   vr(nd,nv)        Data values
!   tmin,tmax        Trimming limits
!   nlag             Number of lags to calculate
!   xlag             Length of the unit lag
!   xltol            Distance tolerance (if <0 then set to xlag/2)
!   ndir             Number of directions to consider
!   azm(ndir)        Azimuth angle of direction (measured positive
!                      degrees clockwise from NS).
!   atol(ndir)       Azimuth (half window) tolerances
!   bandwh           Maximum Horizontal bandwidth (i.e., the deviation 
!                      perpendicular to the defined azimuth).
!   dip(ndir)        Dip angle of direction (measured in negative
!                      degrees down from horizontal).
!   dtol(ndir)       Dip (half window) tolerances
!   bandwd           Maximum "Vertical" bandwidth (i.e., the deviation
!                      perpendicular to the defined dip).
!   isill            1=attempt to standardize, 0=do not
!   sills            the sills (variances) to standardize with
!   nvarg            Number of variograms to compute
!   ivtail(nvarg)    Variable for the tail of each variogram
!   ivhead(nvarg)    Variable for the head of each variogram
!   ivtype(nvarg)    Type of variogram to compute:
!                      1. semivariogram
!                      2. cross-semivariogram
!                      3. covariance
!                      4. correlogram
!                      5. general relative semivariogram
!                      6. pairwise relative semivariogram
!                      7. semivariogram of logarithms
!                      8. rodogram
!                      9. indicator semivariogram (continuous)
!                     10. indicator semivariogram (categorical)
!
!
!
! OUTPUT VARIABLES:  The following arrays are ordered by direction,
!                    then variogram, and finally lag, i.e.,
!                      iloc = (id-1)*nvarg*MAXLG+(iv-1)*MAXLG+il
!
!   np()             Number of pairs
!   dis()            Distance of pairs falling into this lag
!   gam()            Semivariogram, covariance, correlogram,... value
!   hm()             Mean of the tail data
!   tm()             Mean of the head data
!   hv()             Variance of the tail data
!   tv()             Variance of the head data
!
!
!
! PROGRAM NOTES:
!
!   1. The file "gamv.inc" contains the dimensioning parameters.
!      These may have to be changed depending on the amount of data
!      and the requirements to compute different variograms.
!
!
!
! Original:  A.G. Journel                                           1978
! Revisions: K. Guertin                                             1980
!-----------------------------------------------------------------------
      use geostat_gamv
      parameter(PI=3.14159265)
      real      uvxazm(100),uvyazm(100),uvzdec(100), &
                uvhdec(100),csatol(100),csdtol(100)
      logical   omni
!
! Define the distance tolerance if it isn't already:
!
      if(xltol.le.0.0) xltol = 0.5 * xlag
!
! Define the angles and tolerance for each direction:
!
      do id=1,ndir
!
! The mathematical azimuth is measured counterclockwise from EW and
! not clockwise from NS as the conventional azimuth is:
!
            azmuth     = (90.0-azm(id))*PI/180.0
            uvxazm(id) = cos(azmuth)
            uvyazm(id) = sin(azmuth)
            if(atol(id).le.0.0) then
                  csatol(id) = cos(45.0*PI/180.0)
            else
                  csatol(id) = cos(atol(id)*PI/180.0)
            endif
!
! The declination is measured positive down from vertical (up) rather
! than negative down from horizontal:
!
            declin     = (90.0-dip(id))*PI/180.0
            uvzdec(id) = cos(declin)
            uvhdec(id) = sin(declin)
            if(dtol(id).le.0.0) then
                  csdtol(id) = cos(45.0*PI/180.0)
            else
                  csdtol(id) = cos(dtol(id)*PI/180.0)
            endif
      end do
!
! Initialize the arrays for each direction, variogram, and lag:
!
      nsiz = ndir*nvarg*(nlag+2)
      do i=1,nsiz
            np(i)  = 0.
            dis(i) = 0.0
            gam(i) = 0.0
            hm(i)  = 0.0
            tm(i)  = 0.0
            hv(i)  = 0.0
            tv(i)  = 0.0
      end do
      dismxs = ((real(nlag) + 0.5 - EPSLON) * xlag) ** 2
!
! MAIN LOOP OVER ALL PAIRS:
!
      irepo = max(1,min((nd/10),1000))
      do 3 i=1,nd
	  if((int(i/irepo)*irepo).eq.i) write(*,103) i,nd
 103  format('   currently on seed point ',i9,' of ',i9)
      do 4 j=i,nd
!
! Definition of the lag corresponding to the current pair:
!
            dx  = x(j) - x(i)
            dy  = y(j) - y(i)
            dz  = z(j) - z(i)
            dxs = dx*dx
            dys = dy*dy
            dzs = dz*dz
            hs  = dxs + dys + dzs
            if(hs.gt.dismxs) go to 4
            if(hs.lt.0.0) hs = 0.0
            h   = sqrt(hs)
!
! Determine which lag this is and skip if outside the defined distance
! tolerance:
!
            if(h.le.EPSLON) then
                  lagbeg = 1
                  lagend = 1
            else
                  lagbeg = -1
                  lagend = -1
                  do ilag=2,nlag+2
                        if(h.ge.(xlag*real(ilag-2)-xltol).and. &
                           h.le.(xlag*real(ilag-2)+xltol)) then
                              if(lagbeg.lt.0) lagbeg = ilag 
                              lagend = ilag 
                        end if
                  end do
                  if(lagend.lt.0) go to 4
            endif
!
! Definition of the direction corresponding to the current pair. All
! directions are considered (overlapping of direction tolerance cones
! is allowed):
!
            do 5 id=1,ndir
!
! Check for an acceptable azimuth angle:
!
                  dxy = sqrt(max((dxs+dys),0.0))
                  if(dxy.lt.EPSLON) then
                        dcazm = 1.0
                  else
                        dcazm = (dx*uvxazm(id)+dy*uvyazm(id))/dxy
                  endif
                  if(abs(dcazm).lt.csatol(id)) go to 5
!
! Check the horizontal bandwidth criteria (maximum deviation 
! perpendicular to the specified direction azimuth):
!
                  band = uvxazm(id)*dy - uvyazm(id)*dx
                  if(abs(band).gt.bandwh(id)) go to 5
!
! Check for an acceptable dip angle:
!
                  if(dcazm.lt.0.0) dxy = -dxy
                  if(lagbeg.eq.1) then
                        dcdec = 0.0
                  else
                        dcdec = (dxy*uvhdec(id)+dz*uvzdec(id))/h
                        if(abs(dcdec).lt.csdtol(id)) go to 5
                  endif
!
! Check the vertical bandwidth criteria (maximum deviation perpendicular
! to the specified dip direction):
!
                  band = uvhdec(id)*dz - uvzdec(id)*dxy
                  if(abs(band).gt.bandwd(id)) go to 5
!
! Check whether or not an omni-directional variogram is being computed:
!
                  omni = .false.
                  if(atol(id).ge.90.0) omni = .true.
!
! This direction is acceptable - go ahead and compute all variograms:
!
                  do 6 iv=1,nvarg
!
! For this variogram, sort out which is the tail and the head value:
!
                      it = ivtype(iv)
                      if(dcazm.ge.0.0.and.dcdec.ge.0.0) then
                            ii = ivtail(iv)
                            vrh   = vr(i,ii)
                            ii = ivhead(iv)
                            vrt   = vr(j,ii)
                            if(omni.or.it.eq.2) then
                                  ii    = ivhead(iv)
                                  vrtpr = vr(i,ii)
                                  ii    = ivtail(iv)
                                  vrhpr = vr(j,ii)
                            endif
                      else
                            ii = ivtail(iv)
                            vrh   = vr(j,ii)
                            ii = ivhead(iv)
                            vrt   = vr(i,ii)
                            if(omni.or.it.eq.2) then
                                  ii    = ivhead(iv)
                                  vrtpr = vr(j,ii)
                                  ii    = ivtail(iv)
                                  vrhpr = vr(i,ii)
                            endif
                      endif
!
! Reject this pair on the basis of missing values:
!
                      if(vrt.lt.tmin.or.vrh.lt.tmin.or. &
                         vrt.gt.tmax.or.vrh.gt.tmax) go to 6
                      if(it.eq.2.and.(vrtpr.lt.tmin.or.vrhpr.lt.tmin.or. &
                                      vrtpr.gt.tmax.or.vrhpr.gt.tmax)) &
                                                     go to 6
!
!             COMPUTE THE APPROPRIATE "VARIOGRAM" MEASURE
!
!
! The Semivariogram:
!
      if(it.eq.1.or.it.eq.5.or.it.ge.9) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               np(ii)  = np(ii)  + 1.
               dis(ii) = dis(ii) + dble(h)
               tm(ii)  = tm(ii)  + dble(vrt)
               hm(ii)  = hm(ii)  + dble(vrh)
               gam(ii) = gam(ii) + dble((vrh-vrt)*(vrh-vrt))
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and. &
                        vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           gam(ii) = gam(ii) + dble((vrhpr-vrtpr)* &
                                                    (vrhpr-vrtpr))
                     endif
               endif
            end do
!
! The Traditional Cross Semivariogram:
!
      else if(it.eq.2) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               np(ii)  = np(ii)  + 1.
               dis(ii) = dis(ii) + dble(h)
               tm(ii)  = tm(ii)  + dble(0.5*(vrt+vrtpr))
               hm(ii)  = hm(ii)  + dble(0.5*(vrh+vrhpr))
               gam(ii) = gam(ii) + dble((vrhpr-vrh)*(vrt-vrtpr))
            end do
!
! The Covariance:
!
      else if(abs(it).eq.3) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               np(ii)  = np(ii)  + 1.
               dis(ii) = dis(ii) + dble(h)
               tm(ii)  = tm(ii)  + dble(vrt)
               hm(ii)  = hm(ii)  + dble(vrh)
               gam(ii) = gam(ii) + dble(vrh*vrt)
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and. &
                        vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           gam(ii) = gam(ii) + dble(vrhpr*vrtpr)
                     endif
               endif
            end do
!
! The Correlogram:
!
      else if(abs(it).eq.4) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               np(ii)  = np(ii)  + 1.
               dis(ii) = dis(ii) + dble(h)
               tm(ii)  = tm(ii)  + dble(vrt)
               hm(ii)  = hm(ii)  + dble(vrh)
               hv(ii)  = hv(ii)  + dble(vrh*vrh)
               tv(ii)  = tv(ii)  + dble(vrt*vrt)
               gam(ii) = gam(ii) + dble(vrh*vrt)
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and. &
                        vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           hv(ii)  = hv(ii)  + dble(vrhpr*vrhpr)
                           tv(ii)  = tv(ii)  + dble(vrtpr*vrtpr)
                           gam(ii) = gam(ii) + dble(vrhpr*vrtpr)
                     endif
               endif
            end do
!
! The Pairwise Relative:
!
      else if(it.eq.6) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               if(abs(vrt+vrh).gt.EPSLON) then
                     np(ii)  = np(ii)  + 1.
                     dis(ii) = dis(ii) + dble(h)
                     tm(ii)  = tm(ii)  + dble(vrt)
                     hm(ii)  = hm(ii)  + dble(vrh)
                     gamma   = 2.0*(vrt-vrh)/(vrt+vrh)
                     gam(ii) = gam(ii) + dble(gamma*gamma)
               endif
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and. &
                        vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                     if(abs(vrtpr+vrhpr).gt.EPSLON) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           gamma   = 2.0*(vrt-vrh)/(vrt+vrh)
                           gam(ii) = gam(ii) + dble(gamma*gamma)
                     endif
                     endif
               endif
            enddo
!
! Variogram of Logarithms:
!
      else if(it.eq.7) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               if(vrt.gt.EPSLON.and.vrh.gt.EPSLON) then
                     np(ii)  = np(ii)  + 1.
                     dis(ii) = dis(ii) + dble(h)
                     tm(ii)  = tm(ii)  + dble(vrt)
                     hm(ii)  = hm(ii)  + dble(vrh)
                     gamma   = alog(vrt)-alog(vrh)
                     gam(ii) = gam(ii) + dble(gamma*gamma)
               endif
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and. &
                        vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                     if(vrtpr.gt.EPSLON.and.vrhpr.gt.EPSLON) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           gamma   = alog(vrt)-alog(vrh)
                           gam(ii) = gam(ii) + dble(gamma*gamma)
                     endif
                     endif
               endif
            end do
!
! Madogram:
!
      else if(it.eq.8) then
            do il=lagbeg,lagend
               ii = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
               np(ii)  = np(ii)  + 1.
               dis(ii) = dis(ii) + dble(h)
               tm(ii)  = tm(ii)  + dble(vrt)
               hm(ii)  = hm(ii)  + dble(vrh)
               gam(ii) = gam(ii) + dble(abs(vrh-vrt))
               if(omni) then
                     if(vrtpr.ge.tmin.and.vrhpr.ge.tmin.and. &
                        vrtpr.lt.tmax.and.vrhpr.lt.tmax) then
                           np(ii)  = np(ii)  + 1.
                           dis(ii) = dis(ii) + dble(h)
                           tm(ii)  = tm(ii)  + dble(vrtpr)
                           hm(ii)  = hm(ii)  + dble(vrhpr)
                           gam(ii) = gam(ii) + dble(abs(vrhpr-vrtpr))
                     endif
               endif
            end do
      endif
!
! Finish loops over variograms, directions, and the double data loops:
!
 6              continue
 5          continue
 4    continue
 3    continue
!
! Get average values for gam, hm, tm, hv, and tv, then compute
! the correct "variogram" measure:
!
      do 7 id=1,ndir
      do 7 iv=1,nvarg
      do 7 il=1,nlag+2
            i = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
            if(np(i).le.0.) go to 7
            rnum   = np(i)
            dis(i) = dis(i) / dble(rnum)
            gam(i) = gam(i) / dble(rnum)
            hm(i)  = hm(i)  / dble(rnum)
            tm(i)  = tm(i)  / dble(rnum)
            hv(i)  = hv(i)  / dble(rnum)
            tv(i)  = tv(i)  / dble(rnum)
            it = ivtype(iv)
!
! Attempt to standardize:
!
            if(isill.eq.1) then
                  if(ivtail(iv).eq.ivhead(iv)) then
                        iii = ivtail(iv)
                        if((it.eq.1.or.it.ge.9).and.sills(iii).gt.0.0) &
                          gam(i) = gam(i) / sills(iii)
                  end if
            end if
!
!  1. report the semivariogram rather than variogram
!  2. report the cross-semivariogram rather than variogram
!  3. the covariance requires "centering"
!  4. the correlogram requires centering and normalizing
!  5. general relative requires division by lag mean
!  6. report the semi(pairwise relative variogram)
!  7. report the semi(log variogram)
!  8. report the semi(madogram)
!
            if(it.eq.1.or.it.eq.2) then
                  gam(i) = 0.5 * gam(i)
            else if(abs(it).eq.3) then
                  gam(i) = gam(i) - hm(i)*tm(i)
                  if(it.lt.0) then
                        if(sills(ivtail(iv)).lt.0.0.or. &
                           sills(ivhead(iv)).lt.0.0) then
                              gam(i) = -999.0
                        else
                              variance = ( sqrt(sills(ivtail(iv))) &
                                       *   sqrt(sills(ivhead(iv))) )
                              gam(i) = variance - gam(i)
                        end if
                  end if
            else if(abs(it).eq.4) then
                  hv(i)  = hv(i)-hm(i)*hm(i)
                  if(hv(i).lt.0.0) hv(i) = 0.0
                  hv(i)  = dsqrt(hv(i))
                  tv(i)  = tv(i)-tm(i)*tm(i)
                  if(tv(i).lt.0.0) tv(i) = 0.0
                  tv(i)  = dsqrt(tv(i))
                  if((hv(i)*tv(i)).lt.EPSLON) then
                        gam(i) = 0.0
                  else
                        gam(i) =(gam(i)-hm(i)*tm(i))/(hv(i)*tv(i))
                  endif
				  if(it.lt.0) gam(i) = 1.0 - gam(i)
!
! Square "hv" and "tv" so that we return the variance:
!
                  hv(i)  = hv(i)*hv(i)
                  tv(i)  = tv(i)*tv(i)
            else if(it.eq.5) then
                  htave  = 0.5*(hm(i)+tm(i))
                  htave  = htave   *   htave
                  if(htave.lt.EPSLON) then
                        gam(i) = 0.0
                  else
                        gam(i) = gam(i)/dble(htave)
                  endif
            else if(it.ge.6) then
                  gam(i) = 0.5 * gam(i)
            endif
 7    continue
      return
      end
 
 
 
      subroutine writeout
!-----------------------------------------------------------------------
!
!                  Write Out the Results of GAMV3
!                  ******************************
!
! An output file will be written which contains each directional
! variogram ordered by direction and then variogram (the directions
! cycle fastest then the variogram number).  For each variogram there
! will be a one line description and then "nlag" lines with:
!
!        a) lag number (increasing from 1 to nlag)
!        b) separation distance
!        c) the "variogram" value
!        d) the number of pairs for the lag
!        e) the mean of the data contributing to the tail
!        f) the mean of the data contributing to the head
!        g) IF the correlogram - variance of tail values
!        h) IF the correlogram - variance of head values
!
!
!
!
!
!-----------------------------------------------------------------------
      use geostat_gamv
      character title*132
      data      lout/1/
!
! Loop over all the variograms that have been computed:
!
      open(lout,file=outfl,status='UNKNOWN')
      do iv=1,nvarg
!
! Construct a title that reflects the variogram type and the variables
! that were used to calculate the variogram:
!
      it = abs(ivtype(iv))
      if(it.eq. 1) title(1:24) = 'Semivariogram           '
      if(it.eq. 2) title(1:24) = 'Cross Semivariogram     '
      if(it.eq. 3) title(1:24) = 'Covariance              '
      if(it.eq. 4) title(1:24) = 'Correlogram             '
      if(it.eq. 5) title(1:24) = 'General Relative        '
      if(it.eq. 6) title(1:24) = 'Pairwise Relative       '
      if(it.eq. 7) title(1:24) = 'Variogram of Logarithms '
      if(it.eq. 8) title(1:24) = 'Semimadogram            '
      if(it.eq. 9) title(1:24) = 'Indicator 1/2 Variogram '
      if(it.eq.10) title(1:24) = 'Indicator 1/2 Variogram '
      write(title(25:62),100) names(ivtail(iv)),names(ivhead(iv))
 100  format('tail:',a12,' head:',a12)
!
! Loop over all the directions (note the direction in the title):
!
      do id=1,ndir
            write(title(62:74),101) id
 101        format('direction ',i2)
            write(lout,'(a74)') title(1:74)
!
! Write out all the lags:
!
            do il=1,nlag+2
                  i = (id-1)*nvarg*(nlag+2)+(iv-1)*(nlag+2)+il
                  nump = int(np(i))
                  if(it.eq.4) then
                        write(lout,102) il,dis(i),gam(i),nump, &
                                        hm(i),tm(i),hv(i),tv(i)
                  else
                        write(lout,102) il,dis(i),gam(i),nump, &
                                        hm(i),tm(i)
                  endif
 102              format(1x,i3,1x,f12.3,1x,f12.5,1x,i8,4(1x,f14.5))
            end do
      end do
!
! End loop over variograms
!
      end do
!
! Finished:
!
      close(lout)
      return
      end
 
 
 
      subroutine check(vrmin,vrmax)
!-----------------------------------------------------------------------
!
!                Error Check and Note Variogram Types
!                ************************************
!
! Go through each variogram type and note the type to the screen and
! report any possible errors.
!
!
!
!
!
!-----------------------------------------------------------------------
      use geostat_gamv
      real      vrmin(*),vrmax(*)
      character title*80
!
! Loop over all the variograms to be computed:
!
      write(*,*)
      do iv=1,nvarg
!
! Note the variogram type and the variables being used:
!
      it = abs(ivtype(iv))
      if(it.eq. 1) title(1:24) = 'Semivariogram          :'
      if(it.eq. 2) title(1:24) = 'Cross Semivariogram    :'
      if(it.eq. 3) title(1:24) = 'Covariance             :'
      if(it.eq. 4) title(1:24) = 'Correlogram            :'
      if(it.eq. 5) title(1:24) = 'General Relative       :'
      if(it.eq. 6) title(1:24) = 'Pairwise Relative      :'
      if(it.eq. 7) title(1:24) = 'Variogram of Logarithms:'
      if(it.eq. 8) title(1:24) = 'Semimadogram           :'
      if(it.eq. 9) title(1:24) = 'Indicator 1/2 Variogram:'
      if(it.eq.10) title(1:24) = 'Indicator 1/2 Variogram:'
      write(title(25:64),100) names(ivtail(iv)),names(ivhead(iv))
 100  format('  tail=',a12,' head=',a12)
      write(*,101) iv,title(1:64)
 101  format(' Variogram ',i2,1x,a64)
!
! Check for possible errors or inconsistencies:
!
      if(it.eq.2) then
            if(ivtail(iv).eq.ivhead(iv)) write(*,201)
 201        format('  WARNING: cross variogram with the same variable!')
      else if(it.eq.5) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,501)
            if(vrmin(ivtail(iv)).lt.0.0.and.vrmax(ivtail(iv)).gt.0.0) &
                  write(*,502)
            if(vrmin(ivhead(iv)).lt.0.0.and.vrmax(ivhead(iv)).gt.0.0) &
                  write(*,502)
 501        format('  WARNING: cross general relative variogram are', &
                   ' difficult to interpret!')
 502        format('  WARNING: there are both positive and negative', &
                   ' values - lag mean could be zero!')
      else if(it.eq.6) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,601)
            if(vrmin(ivtail(iv)).lt.0.0.and.vrmax(ivtail(iv)).gt.0.0) &
                  write(*,602)
            if(vrmin(ivhead(iv)).lt.0.0.and.vrmax(ivhead(iv)).gt.0.0) &
                  write(*,602)
 601        format('  WARNING: cross pairwise relative variogram are', &
                   ' difficult to interpret!')
 602        format('  WARNING: there are both positive and negative', &
                   ' values - pair means could be zero!')
      else if(it.eq.7) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,701)
            if(vrmin(ivtail(iv)).lt.0.0.or.vrmin(ivhead(iv)).lt.0.0) &
            write(*,702)
 701        format('  WARNING: cross logarithmic variograms may be', &
                   ' difficult to interpret!')
 702        format('  WARNING: there are zero or negative', &
                   ' values - logarithm undefined!')
      else if(it.eq.8) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,801)
 801        format('  WARNING: cross rodograms may be difficult to', &
                   ' interpret!')
      else if(it.eq.9) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,901)
 901        format('  WARNING: cross madograms may be difficult to', &
                   ' interpret!')
      endif
!
! END Loop over all variograms:
!
      end do
!
! Finished:
!
      return
      end



      subroutine makepar
!-----------------------------------------------------------------------
!
!                      Write a Parameter File
!                      **********************
!
!
!
!-----------------------------------------------------------------------
      lun = 99
      open(lun,file='gamv.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for GAMV',/, &
             '                  *******************',/,/, &
             'START OF PARAMETERS:')

      write(lun,11)
 11   format('../data/cluster.dat               ', &
             '-file with data')
      write(lun,12)
 12   format('1   2   0                         ', &
             '-   columns for X, Y, Z coordinates')
      write(lun,13)
 13   format('2   3   4                         ', &
             '-   number of variables,col numbers')
      write(lun,14)
 14   format('-1.0e21     1.0e21                ', &
             '-   trimming limits')
      write(lun,15)
 15   format('gamv.out                          ', &
             '-file for variogram output')
      write(lun,16)
 16   format('10                                ', &
             '-number of lags')
      write(lun,17)
 17   format('5.0                               ', &
             '-lag separation distance')
      write(lun,18)
 18   format('3.0                               ', &
             '-lag tolerance')
      write(lun,19)
 19   format('3                                 ', &
             '-number of directions')
      write(lun,20)
 20   format('0.0  90.0 50.0   0.0  90.0  50.0  ', &
             '-azm,atol,bandh,dip,dtol,bandv')
      write(lun,21)
 21   format('0.0  22.5 25.0   0.0  22.5  25.0  ', &
             '-azm,atol,bandh,dip,dtol,bandv')
      write(lun,22)
 22   format('90.  22.5 25.0   0.0  22.5  25.0  ', &
             '-azm,atol,bandh,dip,dtol,bandv')
      write(lun,23)
 23   format('0                                 ', &
             '-standardize sills? (0=no, 1=yes)')
      write(lun,24)
 24   format('3                                 ', &
             '-number of variograms')
      write(lun,25)
 25   format('1   1   1                         ', &
            '-tail var., head var., variogram type')
      write(lun,26)
 26   format('1   2   2                         ', &
             '-tail var., head var., variogram type')
      write(lun,27)
 27   format('2   2   1                         ', &
             '-tail var., head var., variogram type')
      write(lun,40)
 40   format(//,'type 1 = traditional semivariogram',/, &
                '     2 = traditional cross semivariogram',/, &
                '     3 = covariance',/, &
                '     4 = correlogram',/, &
                '     5 = general relative semivariogram',/, &
                '     6 = pairwise relative semivariogram',/, &
                '     7 = semivariogram of logarithms',/, &
                '     8 = semimadogram',/, &
                '     9 = indicator semivariogram - continuous',/, &
                '     10= indicator semivariogram - categorical') 

      close(lun)
      return
      end
