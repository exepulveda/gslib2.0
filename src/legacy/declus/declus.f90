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
!         DECLUS: a three dimensional cell declustering program
!         *****************************************************
!
! See paper in Computers and Geosciences Vol 15 No 3 (1989) pp 325-332
!
! INPUT/OUTPUT Parameters:
!
!   datafl          file with data
!   icolx,y,z,vr    columns for X, Y, Z, and variable
!   tmin,tmax       trimming limits
!   sumfl           file for summary output
!   outfl           file for output with data & weights
!   yanis,zanis     Y and Z cell anisotropy (Ysize=size*Yanis)
!   iminmax         0=look for minimum declustered mean (1=max)
!   ncell,cmin,cmax number of cell sizes, min size, max size
!   noff            number of origin offsets
!
!
!
! PROGRAM NOTES:
!
!   1. The MAXDAT parameter controls the maximum number of data that may
!      be considered at one time.
!
!   2. The MAXCEL parameter controls the maximum number of cells that
!      may be considered for a particular cell size at one time.
!
!   3. This program requires knowledge of whether the samples are
!      clustered in high or low values or, alternately, knowledge
!      of a ``natural'' cell size, e.g., an underlying regular data
!      spacing.
!
!
!
! The following Parameter controls static dimensioning:
!
!   MAXCEL    maximum number of cells.  The number of cells is a
!             function of the cell size and the size of the area of
!             interest.  In many cases a larger minimum cell size will
!             remove the need to increase MAXCEL.
!
!
! AUTHOR: Clayton V. Deutsch                             DATE: 1989-1999
!-----------------------------------------------------------------------
! use       msflib
      parameter(MAXLEN=512,VERSION=2.905)

      real      var(500)
      integer   test
      character datafl*512,sumfl*512,outfl*512,str*512,strlin*512
      logical   min,testfl
      data      xmin/ 1.0e21/,ymin/ 1.0e21/,zmin/ 1.0e21/, &
                xmax/-1.0e21/,ymax/-1.0e21/,zmax/-1.0e21/, &
                lin/1/,lout/2/,lsum/3/
!
! Declaration of dynamic arrays:
!
      real, allocatable    :: x(:),y(:),z(:),wt(:),vr(:),wtopt(:), &
                              cellwt(:)
      integer, allocatable :: index(:)
!
! Note VERSION number before anything else:
!
      write(*,9999) VERSION
 9999 format(/' DECLUS Version: ',f5.3/)
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
      if(str(1:1).eq.' ') str(1:20) = 'declus.par          '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'declus.par          ') then
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
 1    read(lin,'(a4)',end=97) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
!
! Read Input Parameters:
!
      read(lin,'(a512)',err=97) datafl
      call chknam(datafl,512)
      write(*,*) ' data file = ',datafl(1:40)

      read(lin,*,err=97) ix,iy,iz,ivr
      write(*,*) ' columns = ',ix,iy,iz,ivr

      read(lin,*,err=97) tmin,tmax
      write(*,*) ' tmin,tmax = ',tmin,tmax

      read(lin,'(a512)',err=97) sumfl
      call chknam(sumfl,512)
      write(*,*) ' summary file = ',sumfl(1:40)

      read(lin,'(a512)',err=97) outfl 
      call chknam(outfl,512)
      write(*,*) ' output file = ',outfl(1:40)

      read(lin,*,err=97) anisy,anisz
      write(*,*) ' anisotropy = ',anisy,anisz

      read(lin,*,err=97) minmax
      write(*,*) ' minmax flag = ',minmax

      read(lin,*,err=97) ncell,cmin,cmax
      if(ncell.eq.1) cmax = cmin
      write(*,*) ' ncell min max = ',ncell,cmin,cmax

      read(lin,*,err=97) noff
      write(*,*) ' offsets = ',noff

      write(*,*)
      close(lin)
!
! Make sure that we have a data file:
!
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR: ',datafl,' does not exist'
            write(*,*) '       you need a data file! '
            stop
      endif
!
! Open up the input and output files:
!
      open(lin,file=datafl,status='OLD')
      open(lsum,file=sumfl,status='UNKNOWN')
      open(lout,file=outfl,status='UNKNOWN')
!
! Read the header off the data file, find MAXDAT and 
! prepare the output files:
!
      read(lin,*,err=98)
      read(lin,*,err=98) nvari
      do i=1,nvari
            read(lin,*,err=98)
      end do
      maxdat = 0
 20   read(lin,*,end=40,err=98)(var(j),j=1,nvari)
      if(var(ivr).lt.tmin.or.var(ivr).ge.tmax)go to 20
      maxdat = maxdat +1
      go to 20
 40   continue
!
! Allocate the needed memory:
!
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

      allocate (wt(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if

      allocate (vr(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if

      allocate (wtopt(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if

      allocate (index(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ', &
                       'insufficient memory!', test
            stop
      end if

      rewind(lin)
      read(lin,'(a40)',err=98) str
      write(lout,'(a40)') str
      write(lsum,'(a40)') str
      read(lin,*,err=98) nvari
      write(lout,'(i3)') nvari+1
      do i=1,nvari
            read(lin,'(a40)',err=98) str
            write(lout,'(a40)')      str
      end do
      write(lout,101)
 101  format('Declustering Weight')
      write(lsum,102)
 102  format('2',/,'Cell Size',/,'Declustered Mean')
!
! Now, read in the actual data:
!
      nt = 0
      nd = 0
      vrmin =  1.0e21
      vrmax = -1.0e21
 3    read(lin,*,end=4,err=98) (var(j),j=1,nvari)
!
! Trim this data?
!
      if(var(ivr).lt.tmin.or.var(ivr).ge.tmax) then
            nt = nt + 1
            go to 3
      endif
      nd = nd + 1
!
! Exceeded available storage?
!
      vr(nd) = var(ivr)
      if(ix.le.0) then
            x(nd) = 0.0
      else
            x(nd) = var(ix)
      endif
      if(iy.le.0) then
            y(nd) = 0.0
      else
            y(nd) = var(iy)
      endif
      if(iz.le.0) then
            z(nd) = 0.0
      else
            z(nd) = var(iz)
      endif

      if(vr(nd).lt.vrmin) vrmin = vr(nd)
      if(vr(nd).gt.vrmax) vrmax = vr(nd)

      go to 3
 4    close(lin)
      if(nd.le.1) then
            write(*,*) ' ERROR: there is only one datum'
            stop
      endif
!
! Other Initialization:
!
      min  = .true.
      if(minmax.eq.1) min = .false.
      roff = real(noff)
!
! compute min, max, and average:
!
      vrav = 0.0
      do i=1,nd
            wtopt(i) = 1.0
            vrav  = vrav + vr(i)
            if(x(i).lt.xmin) xmin=x(i)
            if(x(i).gt.xmax) xmax=x(i)
            if(y(i).lt.ymin) ymin=y(i)
            if(y(i).gt.ymax) ymax=y(i)
            if(z(i).lt.zmin) zmin=z(i)
            if(z(i).gt.zmax) zmax=z(i)
      end do
      vrav = vrav / real(nd)
!
! Write Some of the Statistics to the screen:
!
      write(*,900) nd,vrav,vrmin,vrmax,(xmax-xmin), &
                   (ymax-ymin),(zmax-zmin)
 900  format(/' There are ',i8,' data with:',/, &
              '   mean value            = ',f12.5,/, &
              '   minimum and maximum   = ',2f12.5,/, &
              '   size of data vol in X = ',f12.5,/, &
              '   size of data vol in Y = ',f12.5,/, &
              '   size of data vol in Z = ',f12.5,/)
!
! initialize the "best" weight values:
!
      vrop = vrav
      best = 0.0
      write(lsum,300)best,vrop
 300  format(f15.3,2x,f15.3)
!
! define a "lower" origin to use for the cell sizes:
!
      xo1 = xmin - 0.01
      yo1 = ymin - 0.01
      zo1 = zmin - 0.01
!
! define the increment for the cell size:
!
      xinc = (cmax-cmin) / real(ncell)
      yinc = anisy * xinc
      zinc = anisz * xinc
!
! loop over "ncell+1" cell sizes in the grid network:
!
      ncellx = int((xmax-(xo1-cmin))/cmin)+1
      ncelly = int((ymax-(yo1-cmin*anisy))/(cmin*anisy))+1
      ncellz = int((zmax-(zo1-cmin*anisz))/(cmin*anisz))+1
      ncellt = real(ncellx*ncelly*ncellz)
      allocate (cellwt(ncellt),stat = test)
      if (test.ne.0) then
            write(*,*)'Error: Allocation of cell size', &
                        ' failed due to ', &
                        'insufficient memory!', test
            write(*,*)'Trying to allocate ',ncellt
          stop
      end if
      xcs =  cmin        - xinc
      ycs = (cmin*anisy) - yinc
      zcs = (cmin*anisz) - zinc
!
! MAIN LOOP over cell sizes:
!
      do lp=1,ncell+1
            xcs = xcs + xinc
            ycs = ycs + yinc
            zcs = zcs + zinc
!
! initialize the weights to zero:
!
            do i=1,nd
                  wt(i) = 0.0
            end do
!
! determine the maximum number of grid cells in the network:
!
            ncellx = int((xmax-(xo1-xcs))/xcs)+1
            ncelly = int((ymax-(yo1-ycs))/ycs)+1
            ncellz = int((zmax-(zo1-zcs))/zcs)+1
            ncellt = real(ncellx*ncelly*ncellz)
!
! loop over all the origin offsets selected:
!
            xfac = amin1((xcs/roff),(0.5*(xmax-xmin)))
            yfac = amin1((ycs/roff),(0.5*(ymax-ymin)))
            zfac = amin1((zcs/roff),(0.5*(zmax-zmin)))
            do kp=1,noff
                  xo = xo1 - (real(kp)-1.0)*xfac
                  yo = yo1 - (real(kp)-1.0)*yfac
                  zo = zo1 - (real(kp)-1.0)*zfac
!
! initialize the cumulative weight indicators:
!
                  do i=1,ncellt
                        cellwt(i) = 0.0
                  end do
!
! determine which cell each datum is in:
!
                  do i=1,nd
                        icellx = int((x(i) - xo)/xcs) + 1
                        icelly = int((y(i) - yo)/ycs) + 1
                        icellz = int((z(i) - zo)/zcs) + 1
                        icell  = icellx + (icelly-1)*ncellx  &
                                        + (icellz-1)*ncelly*ncellx
                        index(i)      = icell
                        cellwt(icell) = cellwt(icell) + 1.0
                  end do
!
! The weight assigned to each datum is inversely proportional to the
! number of data in the cell.  We first need to get the sum of weights
! so that we can normalize the weights to sum to one:
!
                  sumw = 0.0
                  do i=1,nd
                        ipoint = index(i)
                        sumw   = sumw + (1.0 / cellwt(ipoint))
                  end do
                  sumw = 1.0 / sumw
!
! Accumulate the array of weights (that now sum to one):
!
                  do i=1,nd
                        ipoint = index(i)
                        wt(i) = wt(i) + (1.0/cellwt(ipoint))*sumw
                  end do
!
! End loop over all offsets:
!
            end do
!
! compute the weighted average for this cell size:
!
            sumw  = 0.0
            sumwg = 0.0
            do i=1,nd
                  sumw  = sumw  + wt(i)
                  sumwg = sumwg + wt(i)*vr(i)
            end do
            vrcr  = sumwg / sumw
            write(lsum,300)xcs,vrcr
!
! see if this weighting is optimal:
!
            if((min.and.vrcr.lt.vrop).or.(.not.min.and.vrcr.gt.vrop).or. &
               (ncell.eq.1)) then
                  best = xcs
                  vrop = vrcr
                  do i=1,nd
                        wtopt(i) = wt(i)
                  end do
            endif
!
! END MAIN LOOP over all cell sizes:
!
      end do
      close(3)
!
! Get the optimal weights:
!
      sumw = 0.0
      do i=1,nd
            sumw = sumw + wtopt(i)
      end do
      wtmin = 99999.
      wtmax =-99999.
      facto = real(nd) / sumw
      do i = 1,nd
            wtopt(i) = wtopt(i) * facto
            if(wtopt(i).lt.wtmin) wtmin = wtopt(i)
            if(wtopt(i).gt.wtmax) wtmax = wtopt(i)
      end do
!
! Read the header off the data file:
!
      open(lin,file=datafl,status='OLD')
      read(lin,'()',err=98)
      read(lin,*,err=98) nvari
      do i=1,nvari
            read(lin,'()',err=98)
      end do
!
! Read through input file appending the declustering weight to the
!  output file:
!
      nd2 = 0
 5    read(lin,*,end=6) (var(i),i=1,nvari)
      if(var(ivr).lt.tmin.or.var(ivr).ge.tmax) then
            thewt = -999.0
      else
            nd2   = nd2 + 1
            thewt = wtopt(nd2)
      end if
!
! Write out the results:
!
      backspace lin
      read(lin,'(a)') strlin
      call strlen(strlin,MAXLEN,lostr)
      write(lout,'(a,1x,f10.5)') strlin(1:lostr),thewt
      go to 5
 6    continue
      if(nd.ne.nd2) then
            write(*,*)
            write(*,*) 'ERROR in data somewhere - changed during run!'
            write(*,*)
      end if
!
! Some Debugging Information:
!
      write(*,901) vrop,wtmin,wtmax,1.0
 901  format('   declustered mean      = ',f12.5,/, &
             '   min and max weight    = ',2f12.5,/, &
             '   equal weighting       = ',f12.5,/)
!
! Finished:
!
      close(lin)
      close(lout)
      close(lsum)
      write(*,9998) VERSION
 9998 format(/' DECLUS Version: ',f5.3, ' Finished'/)
      stop
 97   stop 'ERROR in parameter file'
 98   stop 'ERROR in data file'
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
      open(lun,file='declus.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for DECLUS',/, &
             '                  *********************',/,/, &
             'START OF PARAMETERS:')

      write(lun,11)
 11   format('../data/cluster.dat         ', &
             '-file with data')
      write(lun,12)
 12   format('1   2   0   3               ', &
             '-  columns for X, Y, Z, and variable')
      write(lun,13)
 13   format('-1.0e21     1.0e21          ', &
             '-  trimming limits')
      write(lun,14)
 14   format('declus.sum                  ', &
             '-file for summary output')
      write(lun,15)
 15   format('declus.out                  ', &
             '-file for output with data & weights')
      write(lun,16)
 16   format('1.0   1.0                   ', &
             '-Y and Z cell anisotropy (Ysize=size*Yanis)')
      write(lun,17)
 17   format('0                           ', &
             '-0=look for minimum declustered mean (1=max)')
      write(lun,18)
 18   format('24  1.0  25.0               ', &
             '-number of cell sizes, min size, max size')
      write(lun,19)
 19   format('5                           ', &
             '-number of origin offsets')

      close(lun)
      return
      end
