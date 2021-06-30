!
! Module to declare dynamic arrays in multiple subroutines:
!
      module geostat_sgsim

      real,allocatable      :: x(:),y(:),z(:),vr(:),wt(:), &
                vrtr(:),vrgtr(:),close(:),sec(:),sim(:),lvm(:), &
                tmp(:),order(:),covtab(:,:,:),cnodex(:), &
                cnodey(:),cnodez(:),cnodev(:),vra(:),vrea(:)
      real*8,allocatable    :: r(:),rr(:),s(:),a(:)
      integer,allocatable   :: nisb(:),icnode(:)
      integer*2,allocatable :: ixnode(:),iynode(:),iznode(:), &
                ixsbtosr(:),iysbtosr(:),izsbtosr(:)

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
!                Sequential Gaussian Simulation
!                ******************************
!
! The program is executed with no command line arguments.  The user
! will be prompted for the name of a parameter file.  The parameter
! file is described in the documentation (see the example sgsim.par)
!
! The output file will be a GEOEAS file containing the simulated values
! The file is ordered by x,y,z, and then simulation (i.e., x cycles
! fastest, then y, then z, then simulation number).  The values will be
! backtransformed to the original data values if a normal scores
! transform was performed.
!
!
!
! AUTHOR: Clayton V. Deutsch                             DATE: 1989-1999
!-----------------------------------------------------------------------
      use       geostat_sgsim
      include  'sgsim.inc'
!
! Input/Output units used:
!
      lin  = 1
      lout = 2
      ldbg = 3
      llvm = 4
!
! Read the parameters and data (transform as required):
!
      call readparm(MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXSBX,MAXSBY,MAXSBZ)
!
! Call sgsim for the simulation(s):
!
      call sgsim   (MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXSBX,MAXSBY,MAXSBZ)
!
! Finished:
!
      write(*,9998) VERSION
 9998 format(/' SGSIM Version: ',f5.3, ' Finished'/)
      stop
      end



      subroutine readparm(MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ,MAXSBX, &
                          MAXSBY,MAXSBZ)
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
!!       use       msflib
      use       geostat_sgsim
      include  'sgsim.inc'
      real      var(50)
      real*8    p,acorni,cp,oldcp,w
      character transfl*512,smthfl*512,tmpfl*512,datafl*512,outfl*512, &
                dbgfl*512,lvmfl*512,str*512
      logical   testfl,testind
!
! Note VERSION number:
!
      write(*,9999) VERSION
 9999 format(/' SGSIM Version: ',f5.3/)
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
      if(str(1:1).eq.' ') str(1:20) = 'sgsim.par           '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'sgsim.par           ') then
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

      read(lin,*,err=98) ixl,iyl,izl,ivrl,iwt,isecvr
      write(*,*) ' input columns = ',ixl,iyl,izl,ivrl,iwt,isecvr

      read(lin,*,err=98) tmin,tmax
      write(*,*) ' trimming limits = ',tmin,tmax

      read(lin,*,err=98) itrans
      write(*,*) ' transformation flag = ',itrans

      read(lin,'(a512)',err=98) transfl
      call chknam(transfl,512)
      write(*,*) ' transformation file = ',transfl(1:40)

      read(lin,*,err=98) ismooth
      write(*,*) ' consider smoothed distribution (1=yes) = ',ismooth

      read(lin,'(a512)',err=98) smthfl
      call chknam(smthfl,512)
      write(*,*) ' file with smoothed distribution = ',smthfl(1:40)

      read(lin,*,err=98) isvr,iswt
      write(*,*) ' columns = ',isvr,iswt

      read(lin,*,err=98) zmin,zmax
      write(*,*) ' data limits (tails) = ',zmin,zmax

      read(lin,*,err=98) ltail,ltpar
      write(*,*) ' lower tail = ',ltail,ltpar

      read(lin,*,err=98) utail,utpar
      write(*,*) ' upper tail = ',utail,utpar

      read(lin,*,err=98) idbg
      write(*,*) ' debugging level = ',idbg

      read(lin,'(a512)',err=98) dbgfl
      call chknam(dbgfl,512)
      write(*,*) ' debugging file = ',dbgfl(1:40)
      open(ldbg,file=dbgfl,status='UNKNOWN')

      read(lin,'(a512)',err=98) outfl
      call chknam(outfl,512)
      write(*,*) ' output file ',outfl(1:40)

      read(lin,*,err=98) nsim
      write(*,*) ' number of realizations = ',nsim

      read(lin,*,err=98) nx,xmn,xsiz
      write(*,*) ' X grid specification = ',nx,xmn,xsiz

      read(lin,*,err=98) ny,ymn,ysiz
      write(*,*) ' Y grid specification = ',ny,ymn,ysiz

      read(lin,*,err=98) nz,zmn,zsiz
      write(*,*) ' Z grid specification = ',nz,zmn,zsiz

      nxy  = nx*ny
      nxyz = nx*ny*nz

      read(lin,*,err=98) ixv(1)
      write(*,*) ' random number seed = ',ixv(1)
      do i=1,1000
             p = acorni(idum)
      end do

      read(lin,*,err=98) ndmin,ndmax
      write(*,*) ' min and max data = ',ndmin,ndmax

      read(lin,*,err=98) nodmax
      write(*,*) ' maximum previous nodes = ',nodmax

      read(lin,*,err=98) sstrat
      write(*,*) ' two-part search flag = ',sstrat
      if(sstrat.eq.1) ndmax = 0

      read(lin,*,err=98) mults,nmult
      write(*,*) ' multiple grid search flag = ',mults,nmult

      read(lin,*,err=98) noct
      write(*,*) ' number of octants = ',noct

      read(lin,*,err=98) radius,radius1,radius2
      write(*,*) ' search radii = ',radius,radius1,radius2
      if(radius.lt.EPSLON) stop 'radius must be greater than zero'
      radsqd = radius  * radius
      sanis1 = radius1 / radius
      sanis2 = radius2 / radius

      read(lin,*,err=98) sang1,sang2,sang3
      write(*,*) ' search anisotropy angles = ',sang1,sang2,sang3

      read(lin,*,err=98) mxctx,mxcty,mxctz
      write(*,*) ' size of covariance lookup = ',mxctx,mxcty,mxctz
      
      read(lin,*,err=98) ktype
      write(*,*) ' kriging type = ',ktype

      colocorr = 0.0
      if(ktype.eq.4) then
            backspace lin
            read(lin,*,err=98) ktype,colocorr
            varred = 1.0
            backspace lin
            read(lin,*,err=9990) i,xx,varred
 9990       continue
            write(*,*) ' correlation coefficient = ',colocorr
            write(*,*) ' secondary variable varred = ',varred
      end if

      read(lin,'(a512)',err=98) lvmfl
      call chknam(lvmfl,512)
      write(*,*) ' secondary model file = ',lvmfl(1:40)

      read(lin,*,err=98) icollvm
      write(*,*) ' column in secondary model file = ',icollvm

      read(lin,*,err=98) nst(1),c0(1)
      sill = c0(1)
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
            sill     = sill + cc(i)
            if(it(i).eq.4) then
                  write(*,*) ' A power model is NOT allowed '
                  write(*,*) ' Choose a different model and re start '
                  stop
            endif
            write(*,*) ' it,cc,ang[1,2,3]; ',it(i),cc(i), &
                         ang1(i),ang2(i),ang3(i)
            write(*,*) ' a1 a2 a3: ',aa(i),aa1,aa2
      end do
      write(*,*)
      close(lin)
!
! Find the needed parameters:
!
      MAXCTX = mxctx
      MAXCTY = mxcty
      MAXCTZ = mxctz
      MAXCXY = MAXCTX * MAXCTY
      MAXXYZ = MAXCTX * MAXCTY * MAXCTZ
      MAXX   = nx
      MAXY   = ny
      MAXZ   = nz
      MXYZ   = MAXX * MAXY * MAXZ
      if(MXYZ.lt.100) MXYZ = 100
      MAXNOD = nodmax
      MAXSAM = ndmax
      MAXKR1 = MAXNOD + MAXSAM + 1
      MAXKR2 = MAXKR1 * MAXKR1
      MAXSBX = 1
      if(nx.gt.1)then
            MAXSBX = int(nx/2)
            if(MAXSBX.gt.50)MAXSBX=50
      end if

      MAXSBY = 1
      if(ny.gt.1)then
            MAXSBY = int(ny/2)
            if(MAXSBY.gt.50)MAXSBY=50
      end if

      MAXSBZ = 1
      if(nz.gt.1)then
            MAXSBZ = int(nz/2)
            if(MAXSBZ.gt.50)MAXSBZ=50
      end if

      MAXSB = MAXSBX*MAXSBY*MAXSBZ
!
! Find MAXDAT:
!
      MAXDAT = 100
      inquire(file=datafl,exist=testfl)
      if(testfl)then
            open(lin,file=datafl,status='UNKNOWN')
            read(lin,*,err=98)
            read(lin,*,err=99) nvari
            do i=1,nvari
                  read(lin,*)
            end do
            MAXDAT = 0
 33         read(lin,*,end=66,err=98)(var(j),j=1,nvari)
            MAXDAT = MAXDAT + 1
            go to 33
 66         continue
            rewind(lin)
            close(lin)
      end if

      MAXTMP = 1
      inquire(file=smthfl,exist=testfl)
      if(testfl)then
            open(lin,file=smthfl,status='UNKNOWN')
            read(lin,*,err=98)
            read(lin,*,err=97) nvari
            do i=1,nvari
                  read(lin,*)
            end do
            MAXTMP = 0
 22         read(lin,*,end=55,err=97)(var(j),j=1,nvari)
            MAXTMP = MAXTMP + 1
            go to 22
 55         continue
            rewind(lin)
            close(lin)
      end if
      if(MAXTMP.gt.MAXDAT)MAXDAT = MAXTMP
!
! Allocate the needed memory:
!
      allocate(x(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 1: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(y(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 2: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(z(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 3: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(vr(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 4: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(wt(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 5: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(vrtr(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 6: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(vrgtr(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 7: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(close(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 8: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(sec(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 9: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(sim(MXYZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 10: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(lvm(MXYZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 11: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(tmp(MAXXYZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 12: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      MAXORD = MXYZ
      if(MXYZ.lt.MAXCXY) MAXORD=MAXCXY
      allocate(order(MAXORD),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 13: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(covtab(MAXCTX,MAXCTY,MAXCTZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 14: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(cnodex(MAXNOD),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 15: Allocation failed', &
                        ' due to insufficient memory.'
                        stop
            end if

      allocate(cnodey(MAXNOD),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 16: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(cnodez(MAXNOD),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 17: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(cnodev(MAXNOD),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 18: Allocation failed', &
                        ' due to insufficient memory.'
                        stop
            end if

      allocate(vra(MAXKR1),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 19: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(vrea(MAXKR1),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 20: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(r(MAXKR1),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 21: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(rr(MAXKR1),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 22: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(s(MAXKR1),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 23: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(a(MAXKR2),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 24: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(nisb(MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 25: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(icnode(MAXNOD),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 26: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(ixnode(MAXXYZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 27: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(iynode(MAXXYZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 28: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(iznode(MAXXYZ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 29: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(ixsbtosr(8*MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 30: Allocation failed', &
                        ' due to insufficient memory.'
                  stop
            end if

      allocate(iysbtosr(8*MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 31: Allocation failed', &
                       ' due to insufficient memory.'
                  stop
            end if

      allocate(izsbtosr(8*MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR 32: Allocation failed', &
                       ' due to insufficient memory.'
                  stop
            end if
!
! Warn the user if the sill is different than 1.0:
!
      if(sill.gt.(1.0+EPSLON).or.sill.lt.(1.0-EPSLON)) then
            write(*,*) 'WARNING the sill of your variogram is not 1.0!'
            write(*,*) '        the sill = ',sill
            write(*,*)
      end if
!
! Perform some quick error checking:
!
      testfl = .false.
      if(nx.gt.MAXX.or.ny.gt.MAXY.or.nz.gt.MAXZ) then
            write(*,*) 'ERROR: available grid size: ',MAXX,MAXY,MAXZ
            write(*,*) '       you have asked for : ',nx,ny,nz
            testfl = .true.
      end if
      if(ltail.ne.1.and.ltail.ne.2) then
            write(*,*) 'ERROR invalid lower tail option ',ltail
            write(*,*) '      only allow 1 or 2 - see manual '
            testfl = .true.
      endif
      if(utail.ne.1.and.utail.ne.2.and.utail.ne.4) then
            write(*,*) 'ERROR invalid upper tail option ',ltail
            write(*,*) '      only allow 1,2 or 4 - see manual '
            testfl = .true.
      endif
      if(utail.eq.4.and.utpar.lt.1.0) then
            write(*,*) 'ERROR invalid power for hyperbolic tail',utpar
            write(*,*) '      must be greater than 1.0!'
            testfl = .true.
      endif
      if(ltail.eq.2.and.ltpar.lt.0.0) then
            write(*,*) 'ERROR invalid power for power model',ltpar
            write(*,*) '      must be greater than 0.0!'
            testfl = .true.
      endif
      if(utail.eq.2.and.utpar.lt.0.0) then
            write(*,*) 'ERROR invalid power for power model',utpar
            write(*,*) '      must be greater than 0.0!'
            testfl = .true.
      endif
      if(testfl) stop
!
! Check to make sure the data file exists:
!
      nd = 0
      av = 0.0
      ss = 0.0
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'WARNING data file ',datafl,' does not exist!'
            write(*,*) '   - Hope your intention was to create an ', &
                             'unconditional simulation'
            write(*,*) '   - Resetting ndmin, ndmax, and itrans  to 0 '
            write(*,*) '   - Resetting sstrat to 1 '
            ndmin  = 0
            ndmax  = 0
            sstrat = 1
      end if
!
! Establish the reference histogram for the simulation (provided that
! we have data, and we are transforming the data):
!
      if(itrans.eq.1) then
            write(*,*) 'Setting up transformation table'
!
! Decide which file to use for establishing the transformation table:
!
            if(ismooth.eq.1) then
                  tmpfl  = smthfl
                  icolvr = isvr
                  icolwt = iswt
            else
                  tmpfl  = datafl
                  icolvr = ivrl
                  icolwt = iwt
            end if
            inquire(file=tmpfl,exist=testfl)
            if(.not.testfl) then
                  write(*,*) 'ERROR: ',tmpfl,' does not exist'
                  write(*,*) '       this file is needed! '
                  stop
            endif
!
! Open up the file with reference distribution:
!
            open(lin,file=tmpfl,status='UNKNOWN')
            read(lin,'(a40)',err=98) str(1:40)
            read(lin,*,err=99) nvari
            do i=1,nvari
                  read(lin,*,err=98)
            end do
!
! Now, read in the actual data:
!
            nt     = 0
            ntr    = 0
            twt    = 0.0
 3          read(lin,*,end=4,err=99) (var(j),j=1,nvari)
!
! Trim this data?
!
            if(var(icolvr).lt.tmin.or.var(icolvr).ge.tmax) then
                  nt = nt + 1
                  go to 3
            endif
            ntr = ntr + 1
!
! Exceeded available storage?
!
            if(icolvr.gt.nvari.or.icolwt.gt.nvari) then
                  write(*,*) ' ERROR: too few columns in ref data '
                  stop
            endif
!
! Keep this data: Assign the data value and coordinate location:
!
            vrtr(ntr) = var(icolvr)
            if(icolwt.le.0) then
                  vrgtr(ntr) = 1.0
            else
                  vrgtr(ntr) = var(icolwt)
            endif
            if(vrgtr(ntr).le.0.0) then
                  ntr = ntr - 1
                  nt  = nt  + 1
                  go to 3
            end if
            twt = twt + vrgtr(ntr)
!
! Go back for another datum:
!
            go to 3
 4          close(lin)
            if(ntr.le.1) then
                  write(*,*) 'ERROR: too few data for transformation'
                  stop
            endif
!
! Write transformation table:
!
            open(lout,file=transfl,status='UNKNOWN')
!
! Sort data by value:
!
            istart = 1
            iend   = ntr
            call sortem(istart,iend,vrtr,1,vrgtr,c,d,e,f,g,h)
!
! Compute the cumulative probabilities and write transformation table
!
            twt   = max(twt,EPSLON)
            oldcp = 0.0
            cp    = 0.0
            do j=istart,iend
                  cp =  cp + dble(vrgtr(j)/twt)
                  w  = (cp + oldcp)*0.5
                  call gauinv(w,vrg,ierr)
                  if(ierr.eq.1) vrg = UNEST
                  write(lout,201) vrtr(j),vrg
 201              format(f12.5,1x,f12.5)
                  oldcp =  cp
!
! Now, reset the weight to the normal scores value:
!
                  vrgtr(j) = vrg
            end do
            close(lout)
      end if
!
! Now, read the data if the file exists:
!
      inquire(file=datafl,exist=testfl)
      if(testfl) then
            write(*,*) 'Reading input data'
            open(lin,file=datafl,status='OLD')
            read(lin,*,err=99)
            read(lin,*,err=99) nvari
            do i=1,nvari
                  read(lin,*,err=99)
            end do
            if(ixl.gt.nvari.or.iyl.gt.nvari.or.izl.gt.nvari.or. &
               ivrl.gt.nvari.or.isecvr.gt.nvari.or.iwt.gt.nvari) then
                  write(*,*) 'ERROR: you have asked for a column number'
                  write(*,*) '       greater than available in file'
                  stop
            end if
!
! Read all the data until the end of the file:
!
            twt = 0.0
            nd  = 0
            nt  = 0
 5          read(lin,*,end=6,err=99) (var(j),j=1,nvari)
            if(var(ivrl).lt.tmin.or.var(ivrl).ge.tmax) then
                  nt = nt + 1
                  go to 5
            end if
            nd = nd + 1
!
! Acceptable data, assign the value, X, Y, Z coordinates, and weight:
!
            vr(nd) = var(ivrl)
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
            if(iwt.le.0) then
                  wt(nd) = 1.0
            else
                  wt(nd) = var(iwt)
            endif
            if(isecvr.le.0) then
                  sec(nd) = UNEST
            else
                  sec(nd) = var(isecvr)
            endif
!
! Normal scores transform?
!
            if(itrans.eq.1) then
                  vrr = vr(nd)
                  call locate(vrtr,ntr,1,ntr,vrr,j)
                  j   = min(max(1,j),(ntr-1))
                  vrg = powint(vrtr(j),vrtr(j+1),vrgtr(j),vrgtr(j+1), &
                               vrr,1.0)
                  if(vrg.lt.vrgtr(1)  ) vrg = vrgtr(1)
                  if(vrg.gt.vrgtr(ntr)) vrg = vrgtr(nd)
                  vr(nd) = vrg
            end if
            twt = twt + wt(nd)
            av  = av  + var(ivrl)*wt(nd)
            ss  = ss  + var(ivrl)*var(ivrl)*wt(nd)
            go to 5
 6          close(lin)
!
! Compute the averages and variances as an error check for the user:
!
            av = av / max(twt,EPSLON)
            ss =(ss / max(twt,EPSLON)) - av * av
            write(ldbg,111) nd,nt,av,ss
            write(*,   111) nd,nt,av,ss
 111  format(/,' Data for SGSIM: Number of acceptable data  = ',i8,/, &
               '                 Number trimmed             = ',i8,/, &
               '                 Weighted Average           = ',f12.4,/, &
               '                 Weighted Variance          = ',f12.4,/)
      endif
!
! Read secondary attribute model if necessary:
!
      if(ktype.ge.2) then
            write(*,*) 'Reading secondary attribute file'
            inquire(file=lvmfl,exist=testfl)
            if(.not.testfl) then
                  write(*,104) lvmfl
 104              format('WARNING secondary attribute file ',a40, &
                   ' does not exist!')
                  stop
            end if
            open(llvm,file=lvmfl,status='OLD')
            read(llvm,*,err=97)
            read(llvm,*,err=97) nvaril
            do i=1,nvaril
                  read(llvm,*,err=97)
            end do
            index = 0
             
            av = 0.0
            ss = 0.0
            do iz=1,nz
                  do iy=1,ny
                        do ix=1,nx
                           index = index + 1
                           read(llvm,*,err=97) (var(j),j=1,nvaril)
                           lvm(index) = var(icollvm)
                           sim(index) = real(index)
!
! Do we to transform the secondary variable for a local mean?
!
                           if(ktype.eq.2.and.itrans.eq.1) then
                                 vrr = lvm(index)
                                 call locate(vrtr,ntr,1,ntr,vrr,j)
                                 j   =min(max(1,j),(ntr-1))
                                 vrg =powint(vrtr(j),vrtr(j+1),vrgtr(j), &
                                             vrgtr(j+1),vrr,1.0)
                                 if(vrg.lt.vrgtr(1)  ) vrg = vrgtr(1)
                                 if(vrg.gt.vrgtr(ntr)) vrg = vrgtr(nd)
                                 lvm(index) = vrg
                           end if
                           av = av + var(icollvm)
                           ss = ss + var(icollvm)*var(icollvm)
                        end do
                  end do
            end do
            av = av / max(real(nxyz),1.0)
            ss =(ss / max(real(nxyz),1.0)) - av * av
            write(ldbg,112) nxyz,av,ss
            write(*,   112) nxyz,av,ss
 112  format(/,' Secondary Data: Number of data             = ',i8,/, &
               '                 Equal Weighted Average     = ',f12.4,/, &
               '                 Equal Weighted Variance    = ',f12.4,/)
!
! Do we need to work with data residuals? (Locally Varying Mean)
!
            if(ktype.eq.2) then
                  do i=1,nd
                        call getindx(nx,xmn,xsiz,x(i),ix,testind)
                        call getindx(ny,ymn,ysiz,y(i),iy,testind)
                        call getindx(nz,zmn,zsiz,z(i),iz,testind)
                        index = ix + (iy-1)*nx + (iz-1)*nxy
                        sec(i) = lvm(index)
!
! Calculation of residual moved to krige subroutine: vr(i)=vr(i)-sec(i)
!
                  end do
            end if
!
! Do we need to get an external drift attribute for the data?
!
            if(ktype.eq.3) then
                  do i=1,nd
                        if(sec(i).eq.UNEST) then
                              call getindx(nx,xmn,xsiz,x(i),ix,testind)
                              call getindx(ny,ymn,ysiz,y(i),iy,testind)
                              call getindx(nz,zmn,zsiz,z(i),iz,testind)
                              index = ix + (iy-1)*nx + (iz-1)*nxy
                              sec(i) = lvm(index)
                        end if
                  end do
            end if
!
! Transform the secondary attribute to normal scores?
!
            if(ktype.eq.4) then
                  write(ldbg,113) varred
 113              format(/,' Transforming Secondary Data with', &
                           ' variance reduction of ',f12.4,/)
                  write(*,*) 'Transforming secondary variable'
                  write(*,*)
                  call sortem(1,nxyz,lvm,1,sim,c,d,e,f,g,h)
                  oldcp = 0.0
                  cp    = 0.0
                  do i=1,nxyz
                        cp =  cp + dble(1.0/real(nxyz))
                        w  = (cp + oldcp)/2.0
                        call gauinv(w,lvm(i),ierr)
                      lvm(i) = lvm(i) * varred
                        oldcp  =  cp
                  end do
                  call sortem(1,nxyz,sim,1,lvm,c,d,e,f,g,h)
            end if
      end if
!
! Open the output file:
!
      open(lout,file=outfl,status='UNKNOWN')
      write(lout,210)
 210  format('SGSIM Realizations')
      write(lout,211) 1,nx,ny,nz
 211  format(4(1x,i4))
      write(lout,212)
 212  format('value')
      return
!
! Error in an Input File Somewhere:
!
 97   stop 'ERROR in secondary data file!'
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'
      end



      subroutine sgsim(MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ, &
                        MAXSBX,MAXSBY,MAXSBZ)
!-----------------------------------------------------------------------
!
!           Conditional Simulation of a 3-D Rectangular Grid
!           ************************************************
!
! This subroutine generates 3-D realizations of a Gaussian process with
! a given autocovariance model, and conditional to input Gaussian data.
! The conditional simulation is achieved by sequential simulation of all
! the nodes visited by a random path.
!
!
!
! PROGRAM NOTES:
!
!  1. The three dimensional anisotropy parameters, i.e., of the search
!     ellipse and variogram ranges are described in section 2.3 of the
!     manual.   The variogram parameters are described in the same place
!
!  2. The original data and previously simulated grid nodes can be
!     searched separately.  There can be a different maximum number of
!     each and a minimum number of original data can be specified
!     to restrict simulation beyond the limits of the data.  The
!     closeness of previously simulated grid nodes is measured according
!     to the variogram structural distance.
!
!
!
! INPUT VARIABLES:
!
!   nd               Number of data (no missing values)
!   x,y,z(nd)        coordinates of the data
!   vr(nd)           gaussian data (normal scores)
!
!   nx,ny,nz         Number of blocks in X,Y, and Z
!   xmn,ymn,zmn      Coordinate at the center of the first Block
!   xsiz,ysiz,zsiz   spacing of the grid nodes (block size)
!
!   nsim             number of simulations
!   ktype            =1, ordinary kriging; =0, simple kriging
!   sim              the current realization
!   idbg             integer debugging level (0=none,2=normal,4=serious)
!   ldbg             unit number for the debugging output
!   lout             unit number for the output
!
!   radius           Maximum search radius
!   sang1            Azimuth angle of the principal search direction
!   sang2            Dip angle of the principal search direction
!   sang3            Third rotation angle of the search ellipse
!   sanis1           Anisotropy for the dip angle
!   sanis2           Anisotropy for the plunge angle
!   ndmin            Minimum number of data required before sim
!   ndmax            Maximum number of samples for simulation
!   noct             Maximum number per octant if an octant search is
!                      desired (if <= 0, then no octant search)
!
!   nodmax           Maximum number of previously simulated grid nodes
!                      to consider in the simulation.  The structural
!                      variogram distance is used to identify close ones
!
!   c0               Nugget constant (isotropic).
!   cc(nst)          Multiplicative factor of each nested structure.
!   aa(nst)          Parameter "a" of each nested structure.
!   it(nst)          Type of nested structures (1=sph,2=exp,3=gau,4=pow)
!   ang1(nst)        Azimuth angle for the principal direction
!   ang2(nst)        Dip angle for the principal direction
!   ang3(nst)        Third rotation angle to rotate the two minor
!                      directions around the principal direction
!   anis1(nst)       Anisotropy (radius in minor direction at 90
!                      degrees from "ang1" divided by the principal
!                      radius in direction "ang1")
!   anis2(nst)       Anisotropy (radius in minor direction at 90 degrees
!                      vertical from "ang1" divided by the principal
!                      radius in direction "ang1")
!
!
! OUTPUT VARIABLES:  Simulated Values are written to "lout"
!
!
!
! EXTERNAL REFERENCES:
!
!   super            Sets up the super block search of original data
!   search           Search for nearby data values
!   ctable           Builds a covariance table and "spiral" search
!   srchnd           Search for nearby simulated grid nodes
!   sqdist           computes anisotropic squared distance
!   sortem           sorts multiple arrays in ascending order (separate)
!   cova3            Calculates the covariance given a variogram model
!   krige            Sets up and solves either the SK or OK system
!   ksol             Linear system solver using Gaussian elimination
!
!
!
! Concepts taken from F. Alabert and E. Isaaks
!
!-----------------------------------------------------------------------
      use       geostat_sgsim
      include  'sgsim.inc'
      real      var(10)
      real*8    p,acorni,cp,oldcp,w
      logical   testind
!
! Set up the rotation/anisotropy matrices that are needed for the
! variogram and search.
!
      write(*,*) 'Setting up rotation matrices for variogram and search'
      do is=1,nst(1)
            call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is), &
                        is,MAXROT,rotmat)
      end do
      isrot = MAXNST + 1
      call setrot(sang1,sang2,sang3,sanis1,sanis2,isrot,MAXROT,rotmat)
!
! Set up the super block search:
!
      if(sstrat.eq.0) then
            write(*,*) 'Setting up super block search strategy'
            nsec = 1
            call setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z, &
                         vr,wt,nsec,sec,sec2,sec3,MAXSBX,MAXSBY,MAXSBZ, &
                         nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup, &
                         nzsup,zmnsup,zsizsup)
            call picksup(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup, &
                         isrot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr, &
                         iysbtosr,izsbtosr)
      end if
!
! Set up the covariance table and the spiral search:
!
      call ctable(MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ)
!
! MAIN LOOP OVER ALL THE SIMULAUTIONS:
!
      do isim=1,nsim
!
! Read in the secondary data distribution for this realization:
!
            if(isim.gt.1.and.ktype.eq.4) then
                  write(*,*)
                  write(*,*) ' Reading next secondary model'
                  index = 0
                  do iz=1,nz
                        do iy=1,ny
                              do ix=1,nx
                                 index = index + 1
                                 read(llvm,*,end=977)(var(j),j=1,nvaril)
                                 lvm(index) = var(icollvm)
                                 sim(index) = real(index)
                              end do
                        end do
                  end do
                  write(*,*) ' Building CDF from  secondary model'
                  call sortem(1,nxyz,lvm,1,sim,c,d,e,f,g,h)
                  oldcp = 0.0
                  cp    = 0.0
                  do i=1,nxyz
                        cp =  cp + dble(1.0/real(nxyz))
                        w  = (cp + oldcp)/2.0
                        call gauinv(w,lvm(i),ierr)
                        lvm(i) = lvm(i) * varred
                        oldcp  =  cp
                  end do
                  write(*,*) ' Restoring order of secondary model'
                  call sortem(1,nxyz,sim,1,lvm,c,d,e,f,g,h)
 977              continue
                  write(*,*)
            end if
!
! Work out a random path for this realization:
!
            do ind=1,nxyz
                  sim(ind)   = real(acorni(idum))
                  order(ind) = ind
            end do
!
! The multiple grid search works with multiples of 4 (yes, that is
! somewhat arbitrary):
!
            if(mults.eq.1) then
                  do imult=1,nmult
                        nnz = max(1,nz/(imult*4))
                        nny = max(1,ny/(imult*4))
                        nnx = max(1,nx/(imult*4))
                        jz  = 1
                        jy  = 1
                        jx  = 1
                        do iz=1,nnz
                           if(nnz.gt.1) jz = iz*imult*4
                           do iy=1,nny
                              if(nny.gt.1) jy = iy*imult*4
                              do ix=1,nnx
                                 if(nnx.gt.1) jx = ix*imult*4
                                 index = jx + (jy-1)*nx + (jz-1)*nxy
                                 sim(index) = sim(index) - imult
                              end do
                           end do
                        end do
                  end do
            end if
            call sortem(1,nxyz,sim,1,order,c,d,e,f,g,h)
!
! Initialize the simulation:
!
            do ind=1,nxyz
                  sim(ind) = UNEST
            end do
            write(*,*)
            write(*,*) 'Working on realization number ',isim
!
! Assign the data to the closest grid node:
!
            TINY = 0.0001
            do id=1,nd
                  call getindx(nx,xmn,xsiz,x(id),ix,testind)
                  call getindx(ny,ymn,ysiz,y(id),iy,testind)
                  call getindx(nz,zmn,zsiz,z(id),iz,testind)
                  ind = ix + (iy-1)*nx + (iz-1)*nxy
                  xx  = xmn + real(ix-1)*xsiz
                  yy  = ymn + real(iy-1)*ysiz
                  zz  = zmn + real(iz-1)*zsiz
                  test = abs(xx-x(id)) + abs(yy-y(id)) + abs(zz-z(id))
!
! Assign this data to the node (unless there is a closer data):
!
                  if(sstrat.eq.1) then
                        if(sim(ind).ge.0.0) then
                              id2 = int(sim(ind)+0.5)
                              test2 = abs(xx-x(id2)) + abs(yy-y(id2)) &
                                                     + abs(zz-z(id2))
                              if(test.le.test2) sim(ind) = real(id)
                              write(ldbg,102) id,id2
                        else
                              sim(ind) = real(id)
                        end if
                  end if
!
! Assign a flag so that this node does not get simulated:
!
                  if(sstrat.eq.0.and.test.le.TINY) sim(ind)=10.0*UNEST
            end do
 102        format(' WARNING data values ',2i5,' are both assigned to ', &
                 /,'         the same node - taking the closest')
!
! Now, enter data values into the simulated grid:
!
            do ind=1,nxyz
                  id = int(sim(ind)+0.5)
                  if(id.gt.0) sim(ind) = vr(id)
            end do
            irepo = max(1,min((nxyz/10),10000))
!
! MAIN LOOP OVER ALL THE NODES:
!
            do in=1,nxyz
                  if((int(in/irepo)*irepo).eq.in) write(*,103) in
 103              format('   currently on node ',i9)
!
! Figure out the location of this point and make sure it has
! not been assigned a value already:
!
                  index = int(order(in)+0.5)
                  if(sim(index).gt.(UNEST+EPSLON).or. &
                     sim(index).lt.(UNEST*2.0)) go to 5
                  iz = int((index-1)/nxy) + 1
                  iy = int((index-(iz-1)*nxy-1)/nx) + 1
                  ix = index - (iz-1)*nxy - (iy-1)*nx
                  xx = xmn + real(ix-1)*xsiz
                  yy = ymn + real(iy-1)*ysiz
                  zz = zmn + real(iz-1)*zsiz
!
! Now, we'll simulate the point ix,iy,iz.  First, get the close data
! and make sure that there are enough to actually simulate a value,
! we'll only keep the closest "ndmax" data, and look for previously
! simulated grid nodes:
!
                  if(sstrat.eq.0) then
                        call srchsupr(xx,yy,zz,radsqd,isrot,MAXROT, &
                                rotmat,nsbtosr,ixsbtosr,iysbtosr, &
                                izsbtosr,noct,nd,x,y,z,wt,nisb,nxsup, &
                                xmnsup,xsizsup,nysup,ymnsup,ysizsup, &
                                nzsup,zmnsup,zsizsup,nclose,close, &
                                infoct)
                        if(nclose.lt.ndmin) go to 5
                        if(nclose.gt.ndmax) nclose = ndmax
                  endif
                  call srchnd(ix,iy,iz)
!
! Calculate the conditional mean and standard deviation.  This will be
! done with kriging if there are data, otherwise, the global mean and
! standard deviation will be used:
!
                  if(ktype.eq.2) then
                        gmean = lvm(index)
                  else
                        gmean = 0.0
                  end if
                  if((nclose+ncnode).lt.1) then
                        cmean  = gmean
                        cstdev = 1.0
                  else
!
! Perform the kriging.  Note that if there are fewer than four data
! then simple kriging is prefered so that the variance of the
! realization does not become artificially inflated:
!
                        lktype = ktype
                        if(ktype.eq.1.and.(nclose+ncnode).lt.4)lktype=0
                        call krige(ix,iy,iz,xx,yy,zz,lktype,gmean, &
                                   cmean,cstdev,MAXCTX,MAXCTY,MAXCTZ)
                  endif
!
! Draw a random number and assign a value to this node:
!
                  p = acorni(idum)
                  call gauinv(p,xp,ierr)
                  sim(index) = xp * cstdev + cmean
                  if(idbg.ge.3) write(ldbg,141) p,sim(index)
 141              format(' random number ',f6.4,' realization ',f7.4)
!
! Quick check for far out results:
!
                  if(abs(cmean).gt.5.0.or.abs(cstdev).gt.5.0.or. &
                     abs(sim(index)).gt.6.0) then
                  write(ldbg,104) ix,iy,iz,cmean,cstdev,sim(index)
  104             format('WARNING: grid node location: ',3i5,/, &
                         '         conditional mean:   ',f12.5,/, &
                         '         conditional stdev:  ',f12.5,/, &
                         '         simulated value:    ',f12.5)
                  endif
!
! END MAIN LOOP OVER NODES:
!
 5                continue
            end do
!
! Do we need to reassign the data to the grid nodes?
!
            if(sstrat.eq.0) then
                  do id=1,nd
                        call getindx(nx,xmn,xsiz,x(id),ix,testind)
                        call getindx(ny,ymn,ysiz,y(id),iy,testind)
                        call getindx(nz,zmn,zsiz,z(id),iz,testind)
                        xx  = xmn + real(ix-1)*xsiz
                        yy  = ymn + real(iy-1)*ysiz
                        zz  = zmn + real(iz-1)*zsiz
                        ind = ix + (iy-1)*nx + (iz-1)*nxy
                        test=abs(xx-x(id))+abs(yy-y(id))+abs(zz-z(id))
                        if(test.le.TINY) sim(ind) = vr(id)
                  end do
            end if
!
! Back transform each value and write results:
!
            ne = 0
            av = 0.0
            ss = 0.0
            do ind=1,nxyz
                  simval = sim(ind)
                  if(simval.gt.-9.0.and.simval.lt.9.0) then
                        ne = ne + 1
                        av = av + simval
                        ss = ss + simval*simval
                  end if
                  if(itrans.eq.1.and.simval.gt.(UNEST+EPSLON)) then
                        simval = backtr(simval,ntr,vrtr,vrgtr,zmin, &
                                        zmax,ltail,ltpar,utail,utpar)
                        if(simval.lt.zmin) simval = zmin
                        if(simval.gt.zmax) simval = zmax
                  end if
                  write(lout,'(f12.4)') simval
            end do
            av = av / max(real(ne),1.0)
            ss =(ss / max(real(ne),1.0)) - av * av
            write(ldbg,111) isim,ne,av,ss
            write(*,   111) isim,ne,av,ss
 111        format(/,' Realization ',i3,': number   = ',i8,/, &
                     '                  mean     = ',f12.4, &
                     ' (close to 0.0?)',/, &
                     '                  variance = ',f12.4, &
                     ' (close to gammabar(V,V)? approx. 1.0)',/)
!
! END MAIN LOOP OVER SIMULATIONS:
!
      end do
!
! Return to the main program:
!
      return
      end



      subroutine ctable(MAXNOD,MAXCXY,MAXCTX,MAXCTY,MAXCTZ)
!-----------------------------------------------------------------------
!
!               Establish the Covariance Look up Table
!               **************************************
!
! The idea is to establish a 3-D network that contains the covariance
! value for a range of grid node offsets that should be at as large
! as twice the search radius in each direction.  The reason it has to
! be twice as large as the search radius is because we want to use it
! to compute the data covariance matrix as well as the data-point
! covariance matrix.
!
! Secondly, we want to establish a search for nearby nodes that
! in order of closeness as defined by the variogram.
!
!
!
! INPUT VARIABLES:
!
!   xsiz,ysiz,zsiz  Definition of the grid being considered
!   MAXCTX,Y,Z      Number of blocks in covariance table
!
!   covariance table parameters
!
!
!
! OUTPUT VARIABLES:  covtab()         Covariance table
!
! EXTERNAL REFERENCES:
!
!   sqdist          Computes 3-D anisotropic squared distance
!   sortem          Sorts multiple arrays in ascending order
!   cova3           Computes the covariance according to a 3-D model
!
!
!
!-----------------------------------------------------------------------
      use       geostat_sgsim
      parameter(TINY=1.0e-10)
      include  'sgsim.inc'
      real*8    hsqd,sqdist
!
! Size of the look-up table:
!
      nctx = min(((MAXCTX-1)/2),(nx-1))
      ncty = min(((MAXCTY-1)/2),(ny-1))
      nctz = min(((MAXCTZ-1)/2),(nz-1))
!
! Debugging output:
!
      write(ldbg,*)
      write(ldbg,*) 'Covariance Look up table and search for previously'
      write(ldbg,*) 'simulated grid nodes.  The maximum range in each '
      write(ldbg,*) 'coordinate direction for covariance look up is:'
      write(ldbg,*) '          X direction: ',nctx*xsiz
      write(ldbg,*) '          Y direction: ',ncty*ysiz
      write(ldbg,*) '          Z direction: ',nctz*zsiz
      write(ldbg,*) 'Node Values are not searched beyond this distance!'
      write(ldbg,*)
!
! NOTE: If dynamically allocating memory, and if there is no shortage
!       it would a good idea to go at least as far as the radius and
!       twice that far if you wanted to be sure that all covariances
!       in the left hand covariance matrix are within the table look-up.
!
! Initialize the covariance subroutine and cbb at the same time:
!
      call cova3(0.0,0.0,0.0,0.0,0.0,0.0,1,nst,MAXNST,c0,it,cc,aa, &
                 1,MAXROT,rotmat,cmax,cbb)
!
! Now, set up the table and keep track of the node offsets that are
! within the search radius:
!
      nlooku = 0
      do i=-nctx,nctx
      xx = i * xsiz
      ic = nctx + 1 + i
      do j=-ncty,ncty
      yy = j * ysiz
      jc = ncty + 1 + j
      do k=-nctz,nctz
      zz = k * zsiz
      kc = nctz + 1 + k
            call cova3(0.0,0.0,0.0,xx,yy,zz,1,nst,MAXNST,c0,it,cc,aa, &
                       1,MAXROT,rotmat,cmax,covtab(ic,jc,kc))
            hsqd = sqdist(0.0,0.0,0.0,xx,yy,zz,isrot, &
                               MAXROT,rotmat)
            if(real(hsqd).le.radsqd) then
                  nlooku         = nlooku + 1
!
! We want to search by closest variogram distance (and use the
! anisotropic Euclidean distance to break ties:
!
                  tmp(nlooku)   = - (covtab(ic,jc,kc)-TINY*real(hsqd))
                  order(nlooku) = real((kc-1)*MAXCXY+(jc-1)*MAXCTX+ic)
            endif
      end do
      end do
      end do
!
! Finished setting up the look-up table, now order the nodes such
! that the closest ones, according to variogram distance, are searched
! first. Note: the "loc" array is used because I didn't want to make
! special allowance for 2 byte integers in the sorting subroutine:
!
      call sortem(1,nlooku,tmp,1,order,c,d,e,f,g,h)
      do il=1,nlooku
            loc = int(order(il))
            iz  = int((loc-1)/MAXCXY) + 1
            iy  = int((loc-(iz-1)*MAXCXY-1)/MAXCTX) + 1
            ix  = loc-(iz-1)*MAXCXY - (iy-1)*MAXCTX
            iznode(il) = int(iz)
            iynode(il) = int(iy)
            ixnode(il) = int(ix)
      end do
      if(nodmax.gt.MAXNOD) then
            write(ldbg,*)
            write(ldbg,*) 'The maximum number of close nodes = ',nodmax
            write(ldbg,*) 'this was reset from your specification due '
            write(ldbg,*) 'to storage limitations.'
            nodmax = MAXNOD
      endif
!
! Debugging output if requested:
!
      if(idbg.lt.2) return
      write(ldbg,*)
      write(ldbg,*) 'There are ',nlooku,' nearby nodes that will be '
      write(ldbg,*) 'checked until enough close data are found.'
      write(ldbg,*)
      if(idbg.lt.14) return
      do i=1,nlooku
            xx = (ixnode(i) - nctx - 1) * xsiz
            yy = (iynode(i) - ncty - 1) * ysiz
            zz = (iznode(i) - nctz - 1) * zsiz
            write(ldbg,100) i,xx,yy,zz
      end do
 100  format('Point ',i3,' at ',3f12.4)
!
! All finished:
!
      return
      end



      subroutine srchnd(ix,iy,iz)
!-----------------------------------------------------------------------
!
!               Search for nearby Simulated Grid nodes
!               **************************************
!
! The idea is to spiral away from the node being simulated and note all
! the nearby nodes that have been simulated.
!
!
!
! INPUT VARIABLES:
!
!   ix,iy,iz        index of the point currently being simulated
!   sim             the realization so far
!   nodmax          the maximum number of nodes that we want
!   nlooku          the number of nodes in the look up table
!   i[x,y,z]node    the relative indices of those nodes.
!   [x,y,z]mn       the origin of the global grid netwrok
!   [x,y,z]siz      the spacing of the grid nodes.
!
!
!
! OUTPUT VARIABLES:
!
!   ncnode          the number of close nodes
!   icnode()        the number in the look up table
!   cnode[x,y,z]()  the location of the nodes
!   cnodev()        the values at the nodes
!
!
!
!-----------------------------------------------------------------------
      use       geostat_sgsim
      include  'sgsim.inc'
      integer   ninoct(8)
!
! Consider all the nearby nodes until enough have been found:
!
      ncnode = 0
      if(noct.gt.0) then
            do i=1,8
                  ninoct(i) = 0
            end do
      end if
      do 2 il=2,nlooku
            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il))-nctx-1)
            j = iy + (int(iynode(il))-ncty-1)
            k = iz + (int(iznode(il))-nctz-1)
            if(i.lt. 1.or.j.lt. 1.or.k.lt. 1) go to 2
            if(i.gt.nx.or.j.gt.ny.or.k.gt.nz) go to 2
            ind = i + (j-1)*nx + (k-1)*nxy
            if(sim(ind).gt.UNEST) then
!
! Check the number of data already taken from this octant:
!
                  if(noct.gt.0) then
                        idx = ix - i
                        idy = iy - j
                        idz = iz - k
                        if(idz.gt.0) then
                              iq = 4
                              if(idx.le.0 .and. idy.gt.0) iq = 1
                              if(idx.gt.0 .and. idy.ge.0) iq = 2
                              if(idx.lt.0 .and. idy.le.0) iq = 3
                        else
                              iq = 8
                              if(idx.le.0 .and. idy.gt.0) iq = 5
                              if(idx.gt.0 .and. idy.ge.0) iq = 6
                              if(idx.lt.0 .and. idy.le.0) iq = 7
                        end if
                        ninoct(iq) = ninoct(iq) + 1
                        if(ninoct(iq).gt.noct) go to 2
                  end if
                  ncnode = ncnode + 1
                  icnode(ncnode) = il
                  cnodex(ncnode) = xmn + real(i-1)*xsiz
                  cnodey(ncnode) = ymn + real(j-1)*ysiz
                  cnodez(ncnode) = zmn + real(k-1)*zsiz
                  cnodev(ncnode) = sim(ind)
            endif
 2    continue
!
! Return to calling program:
!
      return
      end



      subroutine krige(ix,iy,iz,xx,yy,zz,lktype,gmean,cmean,cstdev, &
                        MAXCTX,MAXCTY,MAXCTZ)
!-----------------------------------------------------------------------
!
!            Builds and Solves the SK or OK Kriging System
!            *********************************************
!
! INPUT VARIABLES:
!
!   ix,iy,iz        index of the point currently being simulated
!   xx,yy,zz        location of the point currently being simulated
!
!
!
! OUTPUT VARIABLES:
!
!   cmean           kriged estimate
!   cstdev          kriged standard deviation
!
!
!
! EXTERNAL REFERENCES: ksol   Gaussian elimination system solution
!
!
!
!-----------------------------------------------------------------------
      use      geostat_sgsim
      include 'sgsim.inc'
      logical first
!
! Size of the kriging system:
!
      first = .false.
      na    = nclose + ncnode
      if(lktype.eq.0) neq = na
      if(lktype.eq.1) neq = na + 1
      if(lktype.eq.2) neq = na
      if(lktype.eq.3) neq = na + 2
      if(lktype.eq.4) neq = na + 1
!
! Set up kriging matrices:
!
      in=0
      do j=1,na
!
! Sort out the actual location of point "j"
!
            if(j.le.nclose) then
                  index  = int(close(j))
                  x1     = x(index)
                  y1     = y(index)
                  z1     = z(index)
                  vra(j) = vr(index)
                  vrea(j)= sec(index)
                  if(lktype.eq.2) vra(j) = vra(j) - vrea(j)
            else
!
! It is a previously simulated node (keep index for table look-up):
!
                  index  = j-nclose
                  x1     = cnodex(index)
                  y1     = cnodey(index)
                  z1     = cnodez(index)
                  vra(j) = cnodev(index)
                  ind    = icnode(index)
                  ix1    = ix + (int(ixnode(ind))-nctx-1)
                  iy1    = iy + (int(iynode(ind))-ncty-1)
                  iz1    = iz + (int(iznode(ind))-nctz-1)
                  index  = ix1 + (iy1-1)*nx + (iz1-1)*nxy
                  vrea(j)= lvm(index)
                  if(lktype.eq.2) vra(j) = vra(j) - vrea(j)
            endif
            do i=1,j
!
! Sort out the actual location of point "i"
!
                  if(i.le.nclose) then
                        index  = int(close(i))
                        x2     = x(index)
                        y2     = y(index)
                        z2     = z(index)
                  else
!
! It is a previously simulated node (keep index for table look-up):
!
                        index  = i-nclose
                        x2     = cnodex(index)
                        y2     = cnodey(index)
                        z2     = cnodez(index)
                        ind    = icnode(index)
                        ix2    = ix + (int(ixnode(ind))-nctx-1)
                        iy2    = iy + (int(iynode(ind))-ncty-1)
                        iz2    = iz + (int(iznode(ind))-nctz-1)
                  endif
!
! Now, get the covariance value:
!
                  in = in + 1
!
! Decide whether or not to use the covariance look-up table:
!
                  if(j.le.nclose.or.i.le.nclose) then
                        call cova3(x1,y1,z1,x2,y2,z2,1,nst,MAXNST,c0,it, &
                                   cc,aa,1,MAXROT,rotmat,cmax,cov)
                        a(in) = dble(cov)
                  else
!
! Try to use the covariance look-up (if the distance is in range):
!
                        ii = nctx + 1 + (ix1 - ix2)
                        jj = ncty + 1 + (iy1 - iy2)
                        kk = nctz + 1 + (iz1 - iz2)
                        if(ii.lt.1.or.ii.gt.MAXCTX.or. &
                           jj.lt.1.or.jj.gt.MAXCTY.or. &
                           kk.lt.1.or.kk.gt.MAXCTZ) then
                              call cova3(x1,y1,z1,x2,y2,z2,1,nst,MAXNST, &
                                   c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)
                        else
                              cov = covtab(ii,jj,kk)
                        endif
                        a(in) = dble(cov)
                  endif
            end do
!
! Get the RHS value (possibly with covariance look-up table):
!
            if(j.le.nclose) then
                  call cova3(xx,yy,zz,x1,y1,z1,1,nst,MAXNST,c0,it,cc,aa, &
                             1,MAXROT,rotmat,cmax,cov)
                  r(j) = dble(cov)
            else
!
! Try to use the covariance look-up (if the distance is in range):
!
                  ii = nctx + 1 + (ix - ix1)
                  jj = ncty + 1 + (iy - iy1)
                  kk = nctz + 1 + (iz - iz1)
                  if(ii.lt.1.or.ii.gt.MAXCTX.or. &
                     jj.lt.1.or.jj.gt.MAXCTY.or. &
                     kk.lt.1.or.kk.gt.MAXCTZ) then
                        call cova3(xx,yy,zz,x1,y1,z1,1,nst,MAXNST,c0,it, &
                                   cc,aa,1,MAXROT,rotmat,cmax,cov)
                  else
                        cov = covtab(ii,jj,kk)
                  endif
                  r(j) = dble(cov)
            endif
            rr(j) = r(j)
      end do
!
! Addition of OK constraint:
!
      if(lktype.eq.1.or.lktype.eq.3) then
            do i=1,na
                  in    = in + 1
                  a(in) = 1.0
            end do
            in       = in + 1
            a(in)    = 0.0
            r(na+1)  = 1.0
            rr(na+1) = 1.0
      endif
!
! Addition of the External Drift Constraint:
!
      if(lktype.eq.3) then
            edmin =  999999.
            edmax = -999999.
            do i=1,na
                  in    = in + 1
                  a(in) = vrea(i)
                  if(a(in).lt.edmin) edmin = a(in)
                  if(a(in).gt.edmax) edmax = a(in)
            end do
            in       = in + 1
            a(in)    = 0.0
            in       = in + 1
            a(in)    = 0.0
            ind      = ix + (iy-1)*nx + (iz-1)*nxy
            r(na+2)  = dble(lvm(ind))
            rr(na+2) = r(na+2)
            if((edmax-edmin).lt.EPSLON) neq = neq - 1
      endif
!
! Addition of Collocated Cosimulation Constraint:
!
      if(lktype.eq.4) then
            sfmin =  1.0e21
            sfmax = -1.0e21
            do i=1,na
                  in    = in + 1
                  a(in) = dble(colocorr)*r(i)
                  if(a(in).lt.sfmin) sfmin = a(in)
                  if(a(in).gt.sfmax) sfmax = a(in)
            end do
            in    = in + 1
            a(in) = 1.0
            ii    = na + 1
            r(ii) = dble(colocorr)
            rr(ii)= r(ii)
!           if((sfmax-sfmin).lt.EPSLON) neq = neq - 1
      end if
!
! Write out the kriging Matrix if Seriously Debugging:
!
      if(idbg.ge.3) then
            write(ldbg,100) ix,iy,iz
            is = 1
            do i=1,neq
                  ie = is + i - 1
                  write(ldbg,101) i,r(i),(a(j),j=is,ie)
                  is = is + i
            end do
 100        format(/,'Kriging Matrices for Node: ',3i4,' RHS first')
 101        format('    r(',i2,') =',f7.4,'  a= ',99f7.4)
      endif
!
! Solve the Kriging System:
!
      if(neq.eq.1.and.lktype.ne.3) then
            s(1)  = r(1) / a(1)
            ising = 0
      else
            call ksol(1,neq,1,a,r,s,ising)
      endif
!
! Write a warning if the matrix is singular:
!
      if(ising.ne.0) then
            if(idbg.ge.1) then
                  write(ldbg,*) 'WARNING SGSIM: singular matrix'
                  write(ldbg,*) '               for node',ix,iy,iz
            endif
            cmean  = gmean
            cstdev = 1.0
            return
      endif
!
! Compute the estimate and kriging variance.  Recall that kriging type
!     0 = Simple Kriging:
!     1 = Ordinary Kriging:
!     2 = Locally Varying Mean:
!     3 = External Drift:
!     4 = Collocated Cosimulation:
!
      cmean  = 0.0
      cstdev = cbb
      sumwts = 0.0
      do i=1,na
            cmean  = cmean  + real(s(i))*vra(i)
            cstdev = cstdev - real(s(i)*rr(i))
            sumwts = sumwts + real(s(i))
      end do

      if(lktype.eq.1) cstdev = cstdev - real(s(na+1))

      if(lktype.eq.2) cmean  = cmean + gmean

      if(lktype.eq.4) then
            ind    = ix + (iy-1)*nx + (iz-1)*nxy
            cmean  = cmean  + real(s(na+1))*lvm(ind)
            cstdev = cstdev - real(s(na+1) *rr(na+1))
      end if
!
! Error message if negative variance:
!
      if(cstdev.lt.0.0) then
            write(ldbg,*) 'ERROR: Negative Variance: ',cstdev
            cstdev = 0.0
      endif
      cstdev = sqrt(max(cstdev,0.0))
!
! Write out the kriging Weights if Seriously Debugging:
!
      if(idbg.ge.3) then
            do i=1,na
                  write(ldbg,140) i,vra(i),s(i)
            end do
 140        format(' Data ',i4,' value ',f8.4,' weight ',f8.4)
            if(lktype.eq.4) write(ldbg,141) lvm(ind),s(na+1)
 141        format(' Sec Data  value ',f8.4,' weight ',f8.4)
            write(ldbg,142) gmean,cmean,cstdev
 142        format(' Global mean ',f8.4,' conditional ',f8.4, &
                   ' std dev ',f8.4)
      end if
!
! Finished Here:
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
      open(lun,file='sgsim.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for SGSIM',/, &
             '                  ********************',/,/, &
             'START OF PARAMETERS:')

      write(lun,11)
 11   format('../data/cluster.dat           ', &
             '-file with data')
      write(lun,12)
 12   format('1  2  0  3  5  0              ', &
             '-  columns for X,Y,Z,vr,wt,sec.var.')
      write(lun,13)
 13   format('-1.0       1.0e21             ', &
             '-  trimming limits')
      write(lun,14)
 14   format('1                             ', &
             '-transform the data (0=no, 1=yes)')
      write(lun,15)
 15   format('sgsim.trn                     ', &
             '-  file for output trans table')
      write(lun,16)
 16   format('0                             ', &
             '-  consider ref. dist (0=no, 1=yes)')
      write(lun,17)
 17   format('histsmth.out                  ', &
             '-  file with ref. dist distribution')
      write(lun,18)
 18   format('1  2                          ', &
             '-  columns for vr and wt')
      write(lun,19)
 19   format('0.0    15.0                   ', &
             '-  zmin,zmax(tail extrapolation)')
      write(lun,20)
 20   format('1       0.0                   ', &
             '-  lower tail option, parameter')
      write(lun,21)
 21   format('1      15.0                   ', &
             '-  upper tail option, parameter')
      write(lun,22)
 22   format('1                             ', &
             '-debugging level: 0,1,2,3')
      write(lun,23)
 23   format('sgsim.dbg                     ', &
             '-file for debugging output')
      write(lun,24)
 24   format('sgsim.out                     ', &
             '-file for simulation output')
      write(lun,25)
 25   format('1                             ', &
             '-number of realizations to generate')
      write(lun,26)
 26   format('50    0.5    1.0              ', &
             '-nx,xmn,xsiz')
      write(lun,27)
 27   format('50    0.5    1.0              ', &
             '-ny,ymn,ysiz')
      write(lun,28)
 28   format('1     0.5    1.0              ', &
             '-nz,zmn,zsiz')
      write(lun,29)
 29   format('69069                         ', &
             '-random number seed')
      write(lun,30)
 30   format('0     8                       ', &
             '-min and max original data for sim')
      write(lun,31)
 31   format('12                            ', &
             '-number of simulated nodes to use')
      write(lun,32)
 32   format('1                             ', &
             '-assign data to nodes (0=no, 1=yes)')
      write(lun,33)
 33   format('1     3                       ', &
             '-multiple grid search (0=no, 1=yes),num')
      write(lun,34)
 34   format('0                             ', &
             '-maximum data per octant (0=not used)')
      write(lun,35)
 35   format('10.0  10.0  10.0              ', &
             '-maximum search radii (hmax,hmin,vert)')
      write(lun,36)
 36   format(' 0.0   0.0   0.0              ', &
             '-angles for search ellipsoid')
      write(lun,37)
 37   format('51    51    11                ', &
             '-size of covariance lookup table')
      write(lun,38)
 38   format('0     0.60   1.0              ', &
             '-ktype: 0=SK,1=OK,2=LVM,3=EXDR,4=COLC')
      write(lun,39)
 39   format('../data/ydata.dat             ', &
             '-  file with LVM, EXDR, or COLC variable')
      write(lun,40)
 40   format('4                             ', &
             '-  column for secondary variable')
      write(lun,41)
 41   format('1    0.1                      ', &
             '-nst, nugget effect')
      write(lun,42)
 42   format('1    0.9  0.0   0.0   0.0     ', &
             '-it,cc,ang1,ang2,ang3')
      write(lun,43)
 43   format('         10.0  10.0  10.0     ', &
             '-a_hmax, a_hmin, a_vert')

      close(lun)
      return
      end
