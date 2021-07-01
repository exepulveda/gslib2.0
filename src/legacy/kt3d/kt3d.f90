!
! Module to declare dynamic arrays in multiple subroutines:
!
      module geostat
      
      integer,allocatable :: nisb(:),ixsbtosr(:),iysbtosr(:),izsbtosr(:)
      real,allocatable    :: x(:),y(:),z(:),vr(:),ve(:),dh(:),tmp(:), &
                close(:),xa(:),ya(:),za(:),vra(:),vea(:),xdb(:),ydb(:), &
                zdb(:),cut(:),cdf(:)
      real*8,allocatable  :: r(:),rr(:),s(:),a(:)
      
      end module
!
!
!
      program main
!-----------------------------------------------------------------------
!
!             Kriging (SK,OK,KT) of a 3-D Rectangular Grid
!             ********************************************
!
! The program is executed with no command line arguments.  The user
! will be prompted for the name of a parameter file.  The parameter
! file is described in the documentation (see the example kt3d.par)
! and should contain the following information:
!
!
!
! AUTHOR: Clayton V. Deutsch                             DATE: 1989-1999
!-----------------------------------------------------------------------
      use       geostat
      include  'kt3d.inc'
!
! Read the parameters, the data, and open the output files:
!
      call readparm(MAXDIS,MAXSBX,MAXSBY,MAXSBZ)
!
! Call kt3d to krige the grid:
!
      call kt3d(MAXDIS,MAXSBX,MAXSBY,MAXSBZ)
!
! Finished:
!
      close(ldbg)
      close(lout)
      write(*,9998) VERSION
 9998 format(/' KT3D Version: ',f5.3, ' Finished'/)
      stop
      end
 
 
 
      subroutine readparm(MAXDIS,MAXSBX,MAXSBY,MAXSBZ)
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
!!      use       msflib
      use       geostat
      include  'kt3d.inc'
      parameter(MV=20)
      real      var(MV)
      character datafl*512,jackfl*512,extfl*512,outfl*512,dbgfl*512, &
                str*512,title*80
      logical   testfl
!
! FORTRAN Units:
!
      lin   = 1
      ldbg  = 3
      lout  = 4
      lext  = 7
      ljack = 8
!
! Note VERSION number:
!
      write(*,9999) VERSION
 9999 format(/' KT3D Version: ',f5.3/)
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
      if(str(1:1).eq.' ') str(1:20) = 'kt3d.par            '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'kt3d.par            ') then
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

      read(lin,*,err=98) idhl,ixl,iyl,izl,ivrl,iextv
      write(*,*) ' columns = ',idhl,ixl,iyl,izl,ivrl,iextv

      read(lin,*,err=98) tmin,tmax
      write(*,*) ' trimming limits = ',tmin,tmax

      read(lin,*,err=98) koption
      write(*,*) ' kriging option = ',koption

!
! This is an undocumented feature to have kt3d construct an IK-type
! distribution:
!
      iktype = 0
      if(koption.lt.0) then
            iktype  = 1
            koption = -koption
      end if
      if(iktype.eq.1) then

            read(lin,*,err=98) ncut
            write(*,*) ' number of cutoffs = ',ncut
!
! Find the needed parameter:
!
            MAXCUT = ncut
!
! Allocate the needed memory:
!21
            allocate(cut(MAXCUT),stat = test)
                  if(test.ne.0)then
                        write(*,*)'ERROR: Allocation failed due to', &
                              ' insufficient memory.'
                        stop
                  end if
!22
            allocate(cdf(MAXCUT),stat = test)
                  if(test.ne.0)then
                        write(*,*)'ERROR: Allocation failed due to', &
                              ' insufficient memory.'
                        stop
                  end if
!
            read(lin,*,err=98) (cut(i),i=1,ncut)
            write(*,*) ' cutoffs = ',(cut(i),i=1,ncut)

      end if

      read(lin,'(a512)',err=98) jackfl
      call chknam(jackfl,512)
      write(*,*) ' jackknife data file = ',jackfl(1:40)

      read(lin,*,err=98) ixlj,iylj,izlj,ivrlj,iextvj
      write(*,*) ' columns = ',ixlj,iylj,izlj,ivrlj,iextvj

      read(lin,*,err=98) idbg
      write(*,*) ' debugging level = ',idbg

      read(lin,'(a512)',err=98) dbgfl
      call chknam(dbgfl,512)
      write(*,*) ' debugging file = ',dbgfl(1:40)

      read(lin,'(a512)',err=98) outfl
      call chknam(outfl,512)
      write(*,*) ' output file = ',outfl(1:40)

      read(lin,*,err=98) nx,xmn,xsiz
      write(*,*) ' nx, xmn, xsiz = ',nx,xmn,xsiz

      read(lin,*,err=98) ny,ymn,ysiz
      write(*,*) ' ny, ymn, ysiz = ',ny,ymn,ysiz

      read(lin,*,err=98) nz,zmn,zsiz
      write(*,*) ' nz, zmn, zsiz = ',nz,zmn,zsiz

      read(lin,*,err=98) nxdis,nydis,nzdis
      write(*,*) ' block discretization:',nxdis,nydis,nzdis

      read(lin,*,err=98) ndmin,ndmax
      write(*,*) ' ndmin,ndmax = ',ndmin,ndmax

      read(lin,*,err=98) noct
      write(*,*) ' max per octant = ',noct

      read(lin,*,err=98) radius,radius1,radius2
      write(*,*) ' search radii = ',radius,radius1,radius2
      if(radius.lt.EPSLON) stop 'radius must be greater than zero'
      radsqd = radius  * radius
      sanis1 = radius1 / radius
      sanis2 = radius2 / radius

      read(lin,*,err=98) sang1,sang2,sang3
      write(*,*) ' search anisotropy angles = ',sang1,sang2,sang3

      read(lin,*,err=98) ktype,skmean
      write(*,*) ' ktype, skmean =',ktype,skmean

      read(lin,*,err=98) (idrif(i),i=1,9)
      write(*,*) ' drift terms = ',(idrif(i),i=1,9)

      read(lin,*,err=98) itrend
      write(*,*) ' itrend = ',itrend

      read(lin,'(a40)',err=98) extfl
      call chknam(extfl,40)
      write(*,*) ' external drift file = ',extfl(1:40)

      read(lin,*,err=98) iextve
      write(*,*) ' variable in external drift file = ',iextve

      read(lin,*,err=98) nst(1),c0(1)
      write(*,*) ' nst, c0 = ',nst(1),c0(1)

      if(nst(1).le.0) then
            write(*,9997) nst(1)
 9997       format(' nst must be at least 1, it has been set to ',i4,/, &
                   ' The c or a values can be set to zero')
            stop
      endif

      do i=1,nst(1)
            read(lin,*,err=98) it(i),cc(i),ang1(i),ang2(i),ang3(i)
            read(lin,*,err=98) aa(i),aa1,aa2
            anis1(i) = aa1 / max(aa(i),EPSLON)
            anis2(i) = aa2 / max(aa(i),EPSLON)
            write(*,*) ' it,cc,ang[1,2,3]; ',it(i),cc(i), &
                         ang1(i),ang2(i),ang3(i)
            write(*,*) ' a1 a2 a3: ',aa(i),aa1,aa2
            if(it(i).eq.4) then
                  if(aa(i).lt.0.0) stop ' INVALID power variogram'
                  if(aa(i).gt.2.0) stop ' INVALID power variogram'
            end if
      end do

      close(lin)
!
! Find the needed parameters:
!
      MAXDIS = nxdis*nydis*nzdis
      MAXSAM = ndmax + 1
      MAXEQ = MAXSAM + MAXDT + 2
      MAXSBX = 1
      if(nx.gt.1)then
            MAXSBX = int(nx/2.00)
            if(MAXSBX.gt.50)MAXSBX=50
      end if
!
      MAXSBY = 1
      if(ny.gt.1)then
            MAXSBY = int(ny/2.00)
            if(MAXSBY.gt.50)MAXSBY=50
      end if
!
      MAXSBZ = 1
      if(nz.gt.1)then
            MAXSBZ = int(nz/2.00)
            if(MAXSBZ.gt.50)MAXSBZ=50
      end if
!
      MAXSB = MAXSBX*MAXSBY*MAXSBZ
      MXSXY = 4 * MAXSBX * MAXSBY
      MXSX  = 2 * MAXSBX
!
! Allocate the needed memory:
!1
      allocate(nisb(MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!2
      allocate(ixsbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!3
      allocate(iysbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!4
      allocate(izsbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!13
      allocate(xa(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!14
      allocate(ya(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!15
      allocate(za(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!16
      allocate(vra(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!17
      allocate(vea(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!18
      allocate(xdb(MAXDIS),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!19
      allocate(ydb(MAXDIS),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!20
      allocate(zdb(MAXDIS),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!23
      allocate(r(MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!24
      allocate(rr(MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!25
      allocate(s(MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!26
      allocate(a(MAXEQ * MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!
! Perform some quick error checking:
!
      if(ndmax.gt.MAXSAM) stop 'ndmax is too big - modify .inc file'
      if(ktype.eq.3.and.iextv.le.0) stop 'must have external variable'
      if(ixl.le.0.and.nx.gt.1) write(*,*) ' WARNING: ixl=0 and nx>1 ! '
      if(iyl.le.0.and.ny.gt.1) write(*,*) ' WARNING: iyl=0 and ny>1 ! '
      if(izl.le.0.and.nz.gt.1) write(*,*) ' WARNING: izl=0 and nz>1 ! '
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
      title(1:22) = 'KT3D ESTIMATES WITH: '
      open(lin,file=datafl,status='OLD')
      read(lin,*)
      read(lin,*,err=99)       nvari
      do i=1,nvari
            read(lin,*)
      end do
      MAXDAT = 0
 22   read(lin,*,end=33,err=99) (var(j),j=1,nvari)
      if(var(ivrl).lt.tmin.or.var(ivrl).ge.tmax) go to 22
      MAXDAT = MAXDAT + 1
      go to 22
 33   continue
!
! Allocate the needed memory:
!5
      allocate(x(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!6
      allocate(y(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!7
      allocate(z(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!8
      allocate(vr(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!9
      allocate(ve(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!10
      allocate(dh(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!11
      allocate(tmp(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!12
      allocate(close(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to', &
                        ' insufficient memory.'
                  stop
            end if
!
      rewind(lin)
      read(lin,'(a58)') title(23:80)
      read(lin,*,err=99)       nvari
      nd = 0
      av = 0.0
      ss = 0.0
      do i=1,nvari
            read(lin,'(a40)',err=99) str
      end do
!
! Some tests on column numbers:
!
      if(ixl.gt.nvari.or.iyl.gt.nvari.or.izl.gt.nvari.or.ivrl.gt.nvari) &
            then
            write(*,*) 'There are only ',nvari,' columns in input data'
            write(*,*) '  your specification is out of range'
            stop
      end if
!
! Read all the data until the end of the file:
!
 2    read(lin,*,end=3,err=99) (var(j),j=1,nvari)
      vrt = var(ivrl)
      if(vrt.lt.tmin.or.vrt.ge.tmax) go to 2
      nd = nd + 1
      if(nd.gt.MAXDAT) then
            write(*,*) ' ERROR: Exceeded available memory for data'
            stop
      end if
!
! Establish the location of this datum:
!
      if(idhl.le.0) then
            dh(nd) = -99
      else
            dh(nd) = var(idhl)
      endif
      if(ixl.le.0) then
            x(nd) = xmn
      else
            x(nd) = var(ixl)
      endif
      if(iyl.le.0) then
            y(nd) = ymn
      else
            y(nd) = var(iyl)
      endif
      if(izl.le.0) then
            z(nd) = zmn
      else
            z(nd) = var(izl)
      endif
!
! Establish the external drift variable (if needed):
!
      ve(nd) = 1.0
      if(ktype.eq.3.or.ktype.eq.2) then
            ve(nd) = var(iextv)
            if(ve(nd).lt.tmin.or.ve(nd).ge.tmax) then
                  write(*,*) ' External drift variable must be present', &
                             ' at all data locations!'
                  write(*,*) ' Encountered at data number ',nd
                  stop
            end if
      end if
      vr(nd) = vrt
      av     = av + vrt
      ss     = ss + vrt*vrt
      go to 2
 3    close(lin)
!
! Compute the averages and variances as an error check for the user:
!
      av = av / max(real(nd),1.0)
      ss =(ss / max(real(nd),1.0)) - av * av
      write(*,*) 'Data for KT3D: Variable number ',ivrl
      write(*,*) '  Number   = ',nd
      write(*,*) '  Average  = ',av
      write(*,*) '  Variance = ',ss
      if(nd.lt.1) then
            write(*,*) ' ERROR: there are no data'
            stop
      end if
!
! Open the debugging and output files:
!
      open(ldbg,file=dbgfl,status='UNKNOWN')
      open(lout,file=outfl,status='UNKNOWN')
      write(lout,'(a80)') title

      if(iktype.eq.0.and.koption.eq.0) then
           write(lout,201) 2,nx,ny,nz
           write(lout,102)
 102       format('Estimate',/,'EstimationVariance')
      end if
      if(iktype.eq.0.and.koption.ge.1) then
           write(lout,201) 7
           write(lout,103)
 103       format('X',/,'Y',/,'Z',/,'True',/,'Estimate',/, &
                  'EstimationVariance',/,'Error: est-true')
      end if
 201  format(4(1x,i4))

      if(iktype.eq.1) then
            if(koption.eq.0) then
                  write(lout,201) ncut,nx,ny,nz
            else
                  write(lout,201) ncut+1
            end if
            do i=1,ncut
                  write(lout,104) i,cut(i)
 104              format('Threshold: ',i2,' = ',f12.5)
            end do
            if(koption.eq.1) write(lout,105)
 105        format('true value')
      end if
!
! Open the external drift file if needed and position it at the
! first grid node in the file:
!
      if((ktype.eq.2.or.ktype.eq.3).and.koption.eq.0) then
            inquire(file=extfl,exist=testfl)
            if(.not.testfl) then
                  write(*,*) 'ERROR file ',extfl,' does not exist!'
                  stop
            endif
            open(lext,file=extfl,status='UNKNOWN')
            read(lext,'(a40)',err=97) str
            read(lext,*,err=97)       nvari
            do i=1,nvari
                  read(lext,'(a40)',err=97) str
            end do
            if(idbg.ge.3) write(ldbg,100) iextve
 100        format('A secondary variable is being used.  The gridded ' &
                   'file',/,'must have the same grid specifications ' &
                   'as the grid you are kriging.',/,'The external ' &
                   'drift variable was taken from column ',i2)
      endif
!
! Set up for cross validation:
!
      if(koption.eq.1) then
            jackfl = datafl
            idhlj  = idhl
            ixlj   = ixl
            iylj   = iyl
            izlj   = izl
            ivrlj  = ivrl
            iextvj = iextv
      end if
!
! Open the file with the jackknife data?
!
      if(koption.gt.0) then
            inquire(file=jackfl,exist=testfl)
            if(.not.testfl) then
                  write(*,*) 'ERROR file ',jackfl,' does not exist!'
                  stop
            endif
            open(ljack,file=jackfl,status='OLD')
            read(ljack,*,err=96)
            read(ljack,*,err=96) nvarij
            do i=1,nvarij
                  read(ljack,*,err=96)
            end do
      end if
!
! Finished here:
!
      return
!
! Error in an Input File Somewhere:
!
 96   stop 'ERROR in jackknife file!'
 97   stop 'ERROR in external drift file!'
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'
      end



      subroutine kt3d(MAXDIS,MAXSBX,MAXSBY,MAXSBZ)
!-----------------------------------------------------------------------
!
!                Krige a 3-D Grid of Rectangular Blocks
!                **************************************
!
! This subroutine estimates point or block values of one variable by
! simple, ordinary, or kriging with a trend model.  It is also possible
! to estimate the trend directly.
!
!
!
! PROGRAM NOTES:
!
!   1. The data and parameters are passed in common blocks defined
!      in kt3d.inc.  Local storage is allocated in the subroutine
!      for kriging matrices, i.e.,
!         - xa,ya,za,vra   arrays for data within search neighborhood
!         - a,r,rr,s       kriging arrays
!         - xdb,ydb,zdb    relative position of discretization points
!         - cbb            block covariance
!   2. The kriged value and the kriging variance is written to Fortran
!      unit number "lout".
!
!
!
!
! Original:  A.G. Journel and C. Lemmer                             1981
! Revisions: A.G. Journel and C. Kostov                             1984
!-----------------------------------------------------------------------
      use        geostat
      include   'kt3d.inc'
      real*8     cbb
      real       var(20)
      logical    first,fircon,accept
      data       fircon/.true./
!
! Set up the rotation/anisotropy matrices that are needed for the
! variogram and search.  Also compute the maximum covariance for
! the rescaling factor:
!
      write(*,*) 'Setting up rotation matrices for variogram and search'
      radsqd = radius * radius
      PMX    = 999.0
      covmax = c0(1)
      do is=1,nst(1)
            call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is), &
                        is,MAXROT,rotmat)
            if(it(is).eq.4) then
                  covmax = covmax + PMX 
            else
                  covmax = covmax + cc(is)
            endif
      end do
      isrot = MAXNST + 1
      call setrot(sang1,sang2,sang3,sanis1,sanis2,isrot,MAXROT,rotmat)
!
! Finish computing the rescaling factor and stop if unacceptable:
!
      if(radsqd.lt.1.0) then
            resc = 2.0 * radius / max(covmax,0.0001)
      else
            resc =(4.0 * radsqd)/ max(covmax,0.0001)
      endif
      if(resc.le.0.0) then
            write(*,*) 'ERROR KT3D: The rescaling value is wrong ',resc
            write(*,*) '            Maximum covariance: ',covmax
            write(*,*) '            search radius:      ',radius
            stop
      endif
      resc = 1.0 / resc
!
! Set up for super block searching:
!
      write(*,*) 'Setting up super block search strategy'
      nsec = 2
      call setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z, &
                   vr,tmp,nsec,ve,dh,sec3,MAXSBX,MAXSBY,MAXSBZ,nisb, &
                   nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,nzsup, &
                   zmnsup,zsizsup)
      call picksup(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup, &
                   isrot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr, &
                   iysbtosr,izsbtosr)
!
! Compute the number of drift terms, if an external drift is being
! considered then it is one more drift term, if SK is being considered
! then we will set all the drift terms off and mdt to 0):
!
      mdt = 1
      do i=1,9
            if(ktype.eq.0.or.ktype.eq.2) idrif(i) = 0
            if(idrif(i).lt.0.or.idrif(i).gt.1) then
                  write(*,*) 'ERROR KT3D: invalid drift term',idrif(i)
                  stop
            endif
            mdt = mdt + idrif(i)
      end do
      if(ktype.eq.3) mdt = mdt + 1
      if(ktype.eq.0) mdt = 0
      if(ktype.eq.2) mdt = 0
!
! Set up the discretization points per block.  Figure out how many
! are needed, the spacing, and fill the xdb,ydb, and zdb arrays with
! the offsets relative to the block center (this only gets done once):
!
! In all cases the offsets are relative to the lower left corner.
! This is done for rescaling the drift terms in the kriging matrix.
!
      if(nxdis.lt.1) nxdis = 1
      if(nydis.lt.1) nydis = 1
      if(nzdis.lt.1) nzdis = 1
      ndb = nxdis * nydis * nzdis
      if(ndb.gt.MAXDIS) then
            write(*,*) 'ERROR KT3D: Too many discretization points',ndb
            write(*,*) '            Increase MAXDIS or lower n[xyz]dis'
            stop
      endif
      xdis = xsiz  / max(real(nxdis),1.0)
      ydis = ysiz  / max(real(nydis),1.0)
      zdis = zsiz  / max(real(nzdis),1.0)
      i    = 0
      xloc = -0.5*(xsiz+xdis)
      do ix =1,nxdis
            xloc = xloc + xdis
            yloc = -0.5*(ysiz+ydis)
            do iy=1,nydis
                  yloc = yloc + ydis
                  zloc = -0.5*(zsiz+zdis)
                  do iz=1,nzdis
                        zloc = zloc + zdis
                        i = i+1
                        xdb(i) = xloc + 0.5*xsiz
                        ydb(i) = yloc + 0.5*ysiz
                        zdb(i) = zloc + 0.5*zsiz
                  end do
            end do
      end do
!
! Initialize accumulators:
!
      nk    = 0
      xk    = 0.0
      vk    = 0.0
      xkmae = 0.0
      xkmse = 0.0
!
! Calculate Block Covariance. Check for point kriging.
!
      call cova3(xdb(1),ydb(1),zdb(1),xdb(1),ydb(1),zdb(1),1,nst,MAXNST, &
                 c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)
!
! Set the ``unbias'' variable so that the matrix solution is more stable
!
      unbias = cov
      cbb    = dble(cov)
      if(ndb.gt.1) then
            cbb = 0.0
            do i=1,ndb
               do j=1,ndb
                  call cova3(xdb(i),ydb(i),zdb(i),xdb(j),ydb(j),zdb(j), &
                     1,nst,MAXNST,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)
                  if(i.eq.j) cov = cov - c0(1)
                  cbb = cbb + dble(cov)
               end do
            end do
            cbb = cbb/dble(real(ndb*ndb))
      end if
      if(idbg.gt.1) then
            write(ldbg,*) ' '
            write(ldbg,*) 'Block Covariance: ',cbb
            write(ldbg,*) ' '
      end if
!
! Mean values of the drift functions:
!
      do i=1,9
            bv(i) = 0.0
      end do
      do i=1,ndb
            bv(1) = bv(1) + xdb(i)
            bv(2) = bv(2) + ydb(i)
            bv(3) = bv(3) + zdb(i)
            bv(4) = bv(4) + xdb(i)*xdb(i)
            bv(5) = bv(5) + ydb(i)*ydb(i)
            bv(6) = bv(6) + zdb(i)*zdb(i)
            bv(7) = bv(7) + xdb(i)*ydb(i)
            bv(8) = bv(8) + xdb(i)*zdb(i)
            bv(9) = bv(9) + ydb(i)*zdb(i)
      end do  
      do i=1,9
            bv(i) = (bv(i) / real(ndb)) * resc
      end do  
!
! Report on progress from time to time:
!
      if(koption.eq.0) then
            nxy   = nx*ny
            nxyz  = nx*ny*nz
            nloop = nxyz
            irepo = max(1,min((nxyz/10),10000))
      else
            nloop = 10000000
            irepo = max(1,min((nd/10),10000))
      end if
      ddh = 0.0
      write(*,*)
      write(*,*) 'Working on the kriging '
!
! MAIN LOOP OVER ALL THE BLOCKS IN THE GRID:
!
      do index=1,nloop
      if((int(index/irepo)*irepo).eq.index) write(*,103) index
 103  format('   currently on estimate ',i9)
!
! Where are we making an estimate?
!
      if(koption.eq.0) then
            iz   = int((index-1)/nxy) + 1
            iy   = int((index-(iz-1)*nxy-1)/nx) + 1
            ix   = index - (iz-1)*nxy - (iy-1)*nx
            xloc = xmn + real(ix-1)*xsiz
            yloc = ymn + real(iy-1)*ysiz
            zloc = zmn + real(iz-1)*zsiz
      else
            read(ljack,*,err=96,end=2) (var(i),i=1,nvarij)
            ddh  = 0.0
            xloc = xmn
            yloc = ymn
            zloc = zmn
            true = UNEST
            secj = UNEST
            if(idhlj.gt.0)  ddh    = var(idhlj)
            if(ixlj.gt.0)   xloc   = var(ixlj)
            if(iylj.gt.0)   yloc   = var(iylj)
            if(izlj.gt.0)   zloc   = var(izlj)
            if(ivrlj.gt.0)  true   = var(ivrlj)
            if(iextvj.gt.0) extest = var(iextvj)
      end if

!
! Read in the external drift variable for this grid node if needed:
!
      if(ktype.eq.2.or.ktype.eq.3) then
            if(koption.eq.0) then
                  read(lext,*) (var(i),i=1,iextve)
                  extest = var(iextve)
            end if
            if(extest.lt.tmin.or.extest.ge.tmax) then
                  est  = UNEST
                  estv = UNEST
                  go to 1
            end if
            resce  = covmax / max(extest,0.0001)
      endif
!
! Find the nearest samples:
!
      call srchsupr(xloc,yloc,zloc,radsqd,isrot,MAXROT,rotmat,nsbtosr, &
                    ixsbtosr,iysbtosr,izsbtosr,noct,nd,x,y,z,tmp, &
                    nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup, &
                    nzsup,zmnsup,zsizsup,nclose,close,infoct)
!
! Load the nearest data in xa,ya,za,vra,vea:
!
      na = 0
      do i=1,nclose
            ind    = int(close(i)+0.5)
            accept = .true.
            if(koption.ne.0.and. &
               (abs(x(ind)-xloc)+abs(y(ind)-yloc)+ abs(z(ind)-zloc)) &
                                 .lt.EPSLON) accept = .false.
            if(koption.ne.0.and. &
               (abs(dh(ind)-ddh)).lt.EPSLON) accept = .false.
            if(accept) then
                  if(na.lt.ndmax) then
                        na = na + 1
                        xa(na)  = x(ind) - xloc + 0.5*xsiz
                        ya(na)  = y(ind) - yloc + 0.5*ysiz
                        za(na)  = z(ind) - zloc + 0.5*zsiz
                        vra(na) = vr(ind)
                        vea(na) = ve(ind)
                  end if
            end if
      end do
!
! Test number of samples found:
!
      if(na.lt.ndmin) then
            est  = UNEST
            estv = UNEST
            go to 1
      end if
!
! Test if there are enough samples to estimate all drift terms:
!
      if(na.ge.1.and.na.le.mdt) then
            if(fircon) then
                  write(ldbg,999)
                  fircon = .false.
            end if
            est  = UNEST
            estv = UNEST
            go to 1
      end if
 999  format(' Encountered a location where there were too few data ',/, &
             ' to estimate all of the drift terms but there would be',/, &
             ' enough data for OK or SK.   KT3D currently leaves ',/, &
             ' these locations unestimated.',/, &
             ' This message is only written once - the first time.',/)
!
! There are enough samples - proceed with estimation.
!
      if(na.le.1) then
!
! Handle the situation of only one sample:
!
            call cova3(xa(1),ya(1),za(1),xa(1),ya(1),za(1),1,nst,MAXNST, &
                      c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb1)
!
! Establish Right Hand Side Covariance:
!
            if(ndb.le.1) then
                  call cova3(xa(1),ya(1),za(1),xdb(1),ydb(1),zdb(1),1, &
                       nst,MAXNST,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb)
            else
                  cb  = 0.0
                  do i=1,ndb
                        call cova3(xa(1),ya(1),za(1),xdb(i),ydb(i), &
                                   zdb(i),1,nst,MAXNST,c0,it,cc,aa,1, &
                                   MAXROT,rotmat,cmax,cov)
                        cb = cb + cov
                        dx = xa(1) - xdb(i)
                        dy = ya(1) - ydb(i)
                        dz = za(1) - zdb(i)
                        if((dx*dx+dy*dy+dz*dz).lt.EPSLON) cb=cb-c0(1)
                  end do
                  cb = cb / real(ndb)
            end if
!vd
!
! Early bug - always did OK in presence of one data.
!
!vd
            if(ktype.eq.2) skmean = extest
            if(ktype.eq.0.or.ktype.eq.2) then
                  wt   = cb / cb1
                  est  = wt * vra(1) + (1.0-wt) * skmean
                  estv = real(cbb) - wt*cb
            else
                  est  = vra(1)
                  estv = real(cbb) - 2.0*cb + cb1
            end if
            nk   = nk + 1
            xk   = xk + est
            vk   = vk + est*est
            go to 1
      end if
!
! Go ahead and set up the OK portion of the kriging matrix:
!
      neq = mdt+na
!
! Initialize the main kriging matrix:
!
      first = .false.
      do i=1,neq*neq
            a(i) = 0.0
      end do
!
! Fill in the kriging matrix:
!
      do i=1,na
      do j=i,na
            call cova3(xa(i),ya(i),za(i),xa(j),ya(j),za(j),1,nst,MAXNST, &
                       c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)
            a(neq*(i-1)+j) = dble(cov)
            a(neq*(j-1)+i) = dble(cov)
      end do
      end do
!
! Fill in the OK unbiasedness portion of the matrix (if not doing SK):
!
      if(neq.gt.na) then
            do i=1,na
                  a(neq*(i-1)+na+1) = dble(unbias)
                  a(neq*na+i)       = dble(unbias)
            end do
      endif
!
! Set up the right hand side:
!
      do i=1,na
            if(ndb.le.1) then
                  call cova3(xa(i),ya(i),za(i),xdb(1),ydb(1),zdb(1),1, &
                       nst,MAXNST,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb)
            else
                  cb  = 0.0
                  do j=1,ndb
                        call cova3(xa(i),ya(i),za(i),xdb(j),ydb(j), &
                                   zdb(j),1,nst,MAXNST,c0,it,cc,aa,1, &
                                   MAXROT,rotmat,cmax,cov)
                        cb = cb + cov
                        dx = xa(i) - xdb(j)
                        dy = ya(i) - ydb(j)
                        dz = za(i) - zdb(j)
                        if((dx*dx+dy*dy+dz*dz).lt.EPSLON) cb=cb-c0(1)
                  end do
                  cb = cb / real(ndb)
            end if
            r(i) = dble(cb)
      end do
      if(neq.gt.na) r(na+1) = dble(unbias)
!
! Add the additional unbiasedness constraints:
!
      im = na + 1
!
! First drift term (linear in "x"):
!
      if(idrif(1).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(xa(k)*resc)
                  a(neq*(k-1)+im) = dble(xa(k)*resc)
            end do
            r(im) = dble(bv(1))
      endif
!
! Second drift term (linear in "y"):
!
      if(idrif(2).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(ya(k)*resc)
                  a(neq*(k-1)+im) = dble(ya(k)*resc)
            end do
            r(im) = dble(bv(2))
      endif
!
! Third drift term (linear in "z"):
!
      if(idrif(3).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(za(k)*resc)
                  a(neq*(k-1)+im) = dble(za(k)*resc)
            end do
            r(im) = dble(bv(3))
      endif
!
! Fourth drift term (quadratic in "x"):
!
      if(idrif(4).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(xa(k)*xa(k)*resc)
                  a(neq*(k-1)+im) = dble(xa(k)*xa(k)*resc)
            end do
            r(im) = dble(bv(4))
      endif
!
! Fifth drift term (quadratic in "y"):
!
      if(idrif(5).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(ya(k)*ya(k)*resc)
                  a(neq*(k-1)+im) = dble(ya(k)*ya(k)*resc)
            end do
            r(im) = dble(bv(5))
      endif
!
! Sixth drift term (quadratic in "z"):
!
      if(idrif(6).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(za(k)*za(k)*resc)
                  a(neq*(k-1)+im) = dble(za(k)*za(k)*resc)
            end do
            r(im) = dble(bv(6))
      endif
!
! Seventh drift term (quadratic in "xy"):
!
      if(idrif(7).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(xa(k)*ya(k)*resc)
                  a(neq*(k-1)+im) = dble(xa(k)*ya(k)*resc)
            end do
            r(im) = dble(bv(7))
      endif
!
! Eighth drift term (quadratic in "xz"):
!
      if(idrif(8).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(xa(k)*za(k)*resc)
                  a(neq*(k-1)+im) = dble(xa(k)*za(k)*resc)
            end do
            r(im) = dble(bv(8))
      endif
!
! Ninth drift term (quadratic in "yz"):
!
      if(idrif(9).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(ya(k)*za(k)*resc)
                  a(neq*(k-1)+im) = dble(ya(k)*za(k)*resc)
            end do
            r(im) = dble(bv(9))
      endif
!
! External drift term (specified by external variable):
!
      if(ktype.eq.3) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(vea(k)*resce)
                  a(neq*(k-1)+im) = dble(vea(k)*resce)
            end do
            r(im) = dble(extest*resce)
      endif
!
! Copy the right hand side to compute the kriging variance later:
!
      do k=1,neq
            rr(k) = r(k)
      end do
      kadim = neq * neq
      ksdim = neq
      nrhs  = 1
      nv    = 1
!
! If estimating the trend then reset all the right hand side terms=0.0:
!
      if(itrend.ge.1) then
            do i=1,na
                  r(i)  = 0.0
                  rr(i) = 0.0
            end do
      endif
!
! Write out the kriging Matrix if Seriously Debugging:
!
      if(idbg.eq.3) then
            write(ldbg,*) 'Estimating node index : ',ix,iy,iz
            is = 1 - neq
            do i=1,neq
                  is = 1 + (i-1)*neq
                  ie = is + neq - 1
                  write(ldbg,100) i,r(i),(a(j),j=is,ie)
 100              format('    r(',i2,') =',f7.4,'  a= ',9(10f7.4))
            end do
      endif
!
! Solve the kriging system:
!
      call ktsol(neq,nrhs,nv,a,r,s,ising,maxeq)
!
! Compute the solution:
!
      if(ising.ne.0) then
            if(idbg.ge.3) write(ldbg,*) ' Singular Matrix ',ix,iy,iz
            est  = UNEST
            estv = UNEST
      else
            est  = 0.0
            estv = real(cbb)
            if(ktype.eq.2) skmean = extest
            do j=1,neq
                  estv = estv - real(s(j))*rr(j)
                  if(j.le.na) then
                        if(ktype.eq.0) then
                              est = est + real(s(j))*(vra(j)-skmean)
                        else if(ktype.eq.2) then
                              est = est + real(s(j))*(vra(j)-vea(j))
                        else
                              est = est + real(s(j))*vra(j)
                        endif
                  endif
            end do
            if(ktype.eq.0.or.ktype.eq.2) est = est + skmean
            nk   = nk + 1
            xk   = xk + est
            vk   = vk + est*est
!
! Write the kriging weights and data if debugging level is above 2:
!
            if(idbg.ge.2) then
                  write(ldbg,*) '       '
                  write(ldbg,*) 'BLOCK: ',ix,iy,iz,' at ',xloc,yloc,zloc
                  write(ldbg,*) '       '
                  if(ktype.ne.0) &
                  write(ldbg,*) '  Lagrange : ',s(na+1)*unbias
                  write(ldbg,*) '  BLOCK EST: x,y,z,vr,wt '
                  do i=1,na
                        xa(i) = xa(i) + xloc - 0.5*xsiz
                        ya(i) = ya(i) + yloc - 0.5*ysiz
                        za(i) = za(i) + zloc - 0.5*zsiz
                        write(ldbg,'(5f12.3)') xa(i),ya(i),za(i), &
                                               vra(i),s(i)
                  end do
                  write(ldbg,*) '  estimate, variance  ',est,estv
            endif
      endif
!
! END OF MAIN KRIGING LOOP:
!
 1          continue
            if(iktype.eq.0) then
                  if(koption.eq.0) then
                        write(lout,'(f9.3,1x,f9.3)') est,estv
                  else
                        err = UNEST
                        if(true.ne.UNEST.and.est.ne.UNEST) then
                              err=est-true
                              xkmae = xkmae + abs(err)
                              xkmse = xkmse + err*err
                        end if
                        write(lout,'(7(f12.3,1x))') xloc,yloc,zloc,true, &
                                              est,estv,err
                  end if
            else
!
! Work out the IK-type distribution implicit to this data configuration
! and kriging weights:
!
                  do icut=1,ncut
                        cdf(icut) = -1.0
                  end do
                  wtmin = 1.0
                  do i=1,na
                        if(s(i).lt.wtmin) wtmin = s(i)
                  end do
                  sumwt = 0.0
                  do i=1,na
                        s(i)  = s(i) - wtmin
                        sumwt = sumwt + s(i)
                  end do
                  do i=1,na
                        s(i) = s(i) / max(0.00001,sumwt)
                  end do
                  if(na.gt.1.and.sumwt.gt.0.00001) then
                        do icut=1,ncut
                              cdf(icut) = 0.0
                              do i=1,na
                                    if(vra(i).le.cut(icut)) &
                                    cdf(icut)=cdf(icut)+s(i)
                              end do
                        end do
                  end if
                  if(koption.eq.0) then
                        write(lout,'(30(f8.4))') (cdf(i),i=1,ncut)
                  else
                        write(lout,'(30(f8.4))') (cdf(i),i=1,ncut),true
                  end if
            end if
      end do
 2    continue
      if(koption.gt.0) close(ljack)
!
! Write statistics of kriged values:
!
 
      if(nk.gt.0.and.idbg.gt.0) then
            xk    = xk/real(nk)
            vk    = vk/real(nk) - xk*xk
            xkmae = xkmae/real(nk)
            xkmse = xkmse/real(nk)
            write(ldbg,105) nk,xk,vk
            write(*,   105) nk,xk,vk
 105        format(/,'Estimated   ',i8,' blocks ',/, &
                     '  average   ',f9.4,/,'  variance  ',f9.4,/)
            if(koption.ne.0) then
                  write(*,106) xkmae,xkmse
 106              format(/,'  mean error',f9.4,/,'  mean sqd e',f9.4)
            end if
      endif
!
! All finished the kriging:
!
      return
 96   stop 'ERROR in jackknife file!'
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
      open(lun,file='kt3d.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for KT3D',/, &
             '                  *******************',/,/, &
             'START OF PARAMETERS:')

      write(lun,11)
 11   format('../data/cluster.dat              ', &
             '-file with data')
      write(lun,12)
 12   format('0  1  2  0  3  0                 ', &
             '-   columns for DH,X,Y,Z,var,sec var')
      write(lun,13)
 13   format('-1.0e21   1.0e21                 ', &
             '-   trimming limits')
      write(lun,14)
 14   format('0                                ', &
             '-option: 0=grid, 1=cross, 2=jackknife')
      write(lun,15)
 15   format('xvk.dat                          ', &
             '-file with jackknife data')
      write(lun,16)
 16   format('1   2   0    3    0              ', &
             '-   columns for X,Y,Z,vr and sec var')
      write(lun,17)
 17   format('3                                ', &
             '-debugging level: 0,1,2,3')
      write(lun,18)
 18   format('kt3d.dbg                         ', &
             '-file for debugging output')
      write(lun,19)
 19   format('kt3d.out                         ', &
             '-file for kriged output')
      write(lun,20)
 20   format('50   0.5    1.0                  ', &
             '-nx,xmn,xsiz')
      write(lun,21)
 21   format('50   0.5    1.0                  ', &
             '-ny,ymn,ysiz')
      write(lun,22)
 22   format('1    0.5    1.0                  ', &
             '-nz,zmn,zsiz')
      write(lun,23)
 23   format('1    1      1                    ', &
             '-x,y and z block discretization')
      write(lun,24)
 24   format('4    8                           ', &
             '-min, max data for kriging')
      write(lun,25)
 25   format('0                                ', &
             '-max per octant (0-> not used)')
      write(lun,26)
 26   format('20.0  20.0  20.0                 ', &
             '-maximum search radii')
      write(lun,27)
 27   format(' 0.0   0.0   0.0                 ', &
             '-angles for search ellipsoid')
      write(lun,28)
 28   format('0     2.302                      ', &
             '-0=SK,1=OK,2=non-st SK,3=exdrift')
      write(lun,29)
 29   format('0 0 0 0 0 0 0 0 0                ', &
             '-drift: x,y,z,xx,yy,zz,xy,xz,zy')
      write(lun,30)
 30   format('0                                ', &
             '-0, variable; 1, estimate trend')
      write(lun,31)
 31   format('extdrift.dat                     ', &
             '-gridded file with drift/mean')
      write(lun,32)
 32   format('4                                ', &
             '-  column number in gridded file')
      write(lun,33)
 33   format('1    0.2                         ', &
             '-nst, nugget effect')
      write(lun,34)
 34   format('1    0.8  0.0   0.0   0.0        ', &
             '-it,cc,ang1,ang2,ang3')
      write(lun,35)
 35   format('         10.0  10.0  10.0        ', &
             '-a_hmax, a_hmin, a_vert')

      close(lun)
      return
      end
