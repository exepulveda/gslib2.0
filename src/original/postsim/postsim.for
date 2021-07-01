      program main
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                                                      %
C Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %
C Junior University.  All rights reserved.                             %
C                                                                      %
C The programs in GSLIB are distributed in the hope that they will be  %
C useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
C responsibility to anyone for the consequences of using them or for   %
C whether they serve any particular purpose or work at all, unless he  %
C says so in writing.  Everyone is granted permission to copy, modify  %
C and redistribute the programs in GSLIB, but only under the condition %
C that this notice and the above copyright notice remain intact.       %
C                                                                      %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c-----------------------------------------------------------------------
c
c                Post Process Simulated Realizations
c                ***********************************
c
c Reads in a set of simulated realizations and post processes them
c according to user specifications.
c
c      1. Computes the E-type mean
c      2. Given a Z-cutoff the program will compute the probability of
c         exceeding the cutoff and the mean value above (and below)
c         the cutoff.
c      3. Given a CDF value the program will compute the corresponding
c         Z-percentile.
c      4. symmetric "P" probability interval
c      5. conditional variances
c
c
c INPUT/OUTPUT Parameters:
c
c      datafl         input realizations
c      outfl          output summary
c      tmin           missing value code
c      iout,outpar    output option and parameter
c      nx,ny,nz       the number of nodes in each coordinate direction
c      nsim           the number of realizations in datafl
c
c
c
c AUTHOR: Clayton V. Deutsch                             DATE: 1989-1999
c-----------------------------------------------------------------------
      use       msflib
      parameter(EPSLON=1.e-12,UNEST=-999.,VERSION=2.907)
      character datafl*512,outfl*512,str*512
      real      meana,meanb
      integer   test
      logical   testfl
      data      lin/1/,lout/2/
c
c Declare dynamic arrays:
c
      real,allocatable :: var(:,:,:),cut(:),cdf(:)
c
c Note VERSION number:
c
      write(*,9999) VERSION
 9999 format(/' POSTSIM Version: ',f5.3/)
c
c Get the name of the parameter file - try the default name if no input:
c
      do i=1,512
            str(i:i) = ' '
      end do
      call getarg(1,str)
      if(str(1:1).eq.' ')then
            write(*,*) 'Which parameter file do you want to use?'
            read (*,'(a)') str
      end if
      if(str(1:1).eq.' ') str(1:20) = 'postsim.par         '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'postsim.par         ') then
                  write(*,*) '        creating a blank parameter file'
                  call makepar
                  write(*,*)
            end if
            stop
      endif
      open(lin,file=str,status='OLD')
c
c Find Start of Parameters:
c
 1    read(lin,'(a4)',end=97) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
c
c Read Input Parameters:
c
      read(lin,'(a512)',err=97) datafl
      call chknam(datafl,512)
      write(*,*) ' data file with simulations = ',datafl(1:40)

      read(lin,*,err=97) nsim
      write(*,*) ' number of simulations = ',nsim

      read(lin,*,err=97) tmin
      write(*,*) ' lower trimming limit = ',tmin

      read(lin,*,err=97) nx,ny,nz
      write(*,*) ' grid size (nx,ny,nz) = ',nx,ny,nz
c
c Find the parameters, then allocate the needed memory.
c
      allocate(var(nx,ny,nsim),stat = test)
      if(test.ne.0)then
            write(*,*)'ERROR: Allocation failed due to',
     +                ' insufficient memory.',test
            stop
      end if

      read(lin,'(a512)',err=97) outfl
      call chknam(outfl,512)
      write(*,*) ' output file = ',outfl(1:40)

      read(lin,*,err=97) iout,outpar
      write(*,*) ' output option and parameter:',iout,outpar

      if(iout.eq.3) then
            if(outpar.lt.0.0) stop 'Invalid p-value for iout=3'
            if(outpar.gt.1.0) stop 'Invalid p-value for iout=3'
      end if

      close(lin)
c
c Allocation:
c
      allocate(cdf(nsim),stat = test)
      if(test.ne.0)then
            write(*,*)'ERROR: Allocation failed due to',
     +                ' insufficient memory.',test
            stop
      end if
      allocate(cut(nsim),stat = test)
      if(test.ne.0)then
            write(*,*)'ERROR: Allocation failed due to',
     +                ' insufficient memory.',test
            stop
      end if
c
c Set up cdf once:
c
      cdfinc = 1.0/real(nsim)
      cdf(1) = 0.5*cdfinc
      do i=2,nsim
            cdf(i) = cdf(i-1) + cdfinc
      end do
      if(iout.eq.3) then
            if(outpar.lt.cdf(1)   ) outpar = cdf(1)
            if(outpar.gt.cdf(nsim)) outpar = cdf(nsim)
      end if
c
c Open input file with all of the realizations and output file:
c
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) stop 'ERROR datafl does not exist!'
      open(lin,file=datafl,status='UNKNOWN')
      open(lout,file=outfl,status='UNKNOWN')

      if(iout.eq.1) then
            write(lout,101)
 101        format('E-type mean values')
            write(lout,201) 1,nx,ny,nz
            write(lout,102)
 102        format('mean')
      end if
      if(iout.eq.2) then
            write(lout,103) outpar
 103        format('Probability and mean value > ',f12.4)
            write(lout,201) 3,nx,ny,nz
            write(lout,104)
 104        format('prob > cutoff',/,'mean > cutoff',/,'mean < cutoff')
      end if
      if(iout.eq.3) then
            write(lout,105) outpar
 105        format('Z value corresponding to CDF = ',f7.4)
            write(lout,201) 1,nx,ny,nz
            write(lout,106) outpar
 106        format('value')
      end if
      if(iout.eq.4) then
            write(lout,107) outpar
 107        format('Probability Interval = ',f7.4)
            write(lout,201) 2,nx,ny,nz
            write(lout,108)
 108        format('lower',/,'upper')
      end if
      if(iout.eq.5) then
            write(lout,109)
 109        format('Conditional Variance')
            write(lout,201) 1,nx,ny,nz
            write(lout,110)
 110        format('variance')
      end if
 201  format(4(1x,i4))

c
c If we are getting the Z value for a particular probability we can
c establish the weighting straight away:
c
      if(iout.eq.3) then
            call locate(cdf,nsim,1,nsim,outpar,iii)
            wtiii  = (outpar-cdf(iii))/(cdf(iii+1)-cdf(iii))
            wtiiii = 1.0 - wtiii
      endif
      if(iout.eq.4) then
            outlow = (1.0-outpar)/2.0
            call locate(cdf,nsim,1,nsim,outlow,iii)
            wtiii  = (outlow-cdf(iii))/(cdf(iii+1)-cdf(iii))
            wtiiii = 1.0 - wtiii
            outupp = (1.0+outpar)/2.0
            call locate(cdf,nsim,1,nsim,outupp,jjj)
            wtjjj  = (outupp-cdf(jjj))/(cdf(jjj+1)-cdf(jjj))
            wtjjjj = 1.0 - wtjjj
      endif
c
c MAIN LOOP OVER ALL OF THE Z LEVELS:
c
      do 2 iznow=1,nz
c
c Rewind data file and read in the distributions
c
      rewind(lin)
      read(lin,*)
      read(lin,*) nvari
      do i=1,nvari
            read(lin,*)
      end do
      do is=1,nsim
            do iz=1,nz
                  do iy=1,ny
                        do ix=1,nx
                              if(iz.eq.iznow) then
                                    read(lin,*)var(ix,iy,is)
                              else
                                    read(lin,*)
                              endif
                        end do
                  end do
            end do
      end do
c
c Now, we have the distributions to work with. Go through each node:
c
      do iy=1,ny
      do ix=1,nx
c
c Load the cut array and sort:
c
            do is=1,nsim
                  cut(is) = var(ix,iy,is)
            end do
            call sortem(1,nsim,cut,0,b,c,d,e,f,g,h)
c
c Compute the E-type?
c
            if(iout.eq.1) then
                  if(cut(nsim).lt.tmin) then
                        etype = UNEST
                  else
                        etype = 0.0
                        do is=1,nsim
                              etype = etype + cut(is)
                        end do
                        etype = etype / real(nsim)
                  endif
                  write(lout,'(f9.4)') etype
c
c Compute the probability and mean above cutoff?
c
            else if(iout.eq.2) then
                  if(cut(nsim).lt.tmin) then
                        prob  = UNEST
                        meana = UNEST
                        meanb = UNEST
                  else
                        prob  = 0.0
                        meana = 0.0
                        meanb = 0.0
                        do is=1,nsim
                              if(cut(is).ge.outpar) then
                                    prob  = prob + 1.0
                                    meana = meana + cut(is)
                              else
                                    meanb = meanb + cut(is)
                              endif
                        end do
                        if(prob.eq.0) then
                              meana = UNEST
                        else
                              meana = meana / prob
                        endif
                        if((real(nsim)-prob).eq.0) then
                              meanb = UNEST
                        else
                              meanb = meanb / (real(nsim)-prob)
                        endif
                        prob  = prob / real(nsim)
                  endif
                  write(lout,'(f9.4,1x,f9.4,2x,f9.4)') prob,meana,meanb
c
c Now look for the right Z value?
c
            else if(iout.eq.3) then
                  if(cut(nsim).lt.tmin) then
                        zval = UNEST
                  else
                        zval = wtiii*cut(iii)+wtiiii*cut(iii+1)
                  endif
                  write(lout,'(f9.4)') zval
c
c Now look for the right Z values for probability interval:
c
            else if(iout.eq.4) then
                  if(cut(nsim).lt.tmin) then
                        zlow = UNEST
                  else
                        zlow = wtiii*cut(iii)+wtiiii*cut(iii+1)
                  endif
                  if(cut(nsim).lt.tmin) then
                        zupp = UNEST
                  else
                        zupp = wtjjj*cut(jjj)+wtjjjj*cut(jjj+1)
                  endif
                  write(lout,'(f12.4,1x,f12.4)') zlow,zupp
c
c Conditional variance?
c
            else if(iout.eq.5) then
                  if(cut(nsim).lt.tmin) then
                        cvar  = UNEST
                  else
                        etype = 0.0
                        cvar  = 0.0
                        do is=1,nsim
                              etype = etype + cut(is)
                              cvar  = cvar  + cut(is)*cut(is)
                        end do
                        etype = etype / real(nsim)
                        cvar  = cvar  / real(nsim)-etype*etype
                  endif
                  write(lout,'(f12.4)') cvar
            endif
c
c End loop over this ix, iy location:
c
        end do
        end do
c
c End loop over this level and then loop over all levels:
c
 2      continue
c
c Finished:
c
      write(*,9998) VERSION
 9998 format(/' POSTSIM Version: ',f5.3, ' Finished'/)
      stop
 97   stop 'ERROR in parameter file!'
      end



      subroutine makepar
c-----------------------------------------------------------------------
c
c                      Write a Parameter File
c                      **********************
c
c
c
c-----------------------------------------------------------------------
      lun = 99
      open(lun,file='postsim.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for POSTSIM',/,
     +       '                  **********************',/,/,
     +       'START OF PARAMETERS:')

      write(lun,11)
 11   format('sgsim.out                        ',
     +       '-file with simulated realizations')
      write(lun,12)
 12   format('50                               ',
     +       '-   number of realizations')
      write(lun,13)
 13   format('-0.001   1.0e21                  ',
     +       '-   trimming limits')
      write(lun,14)
 14   format('20   20   1                      ',
     +       '-nx, ny, nz')
      write(lun,15)
 15   format('postsim.out                      ',
     +       '-file for output array(s)')
      write(lun,16)
 16   format('2   0.25                         ',
     +       '-output option, output parameter')
      write(lun,17)
 17   format(//,'option 1 = E-type mean',/,
     +          '       2 = prob and mean above threshold (par)',/,
     +          '       3 = Z-percentile corresponding to (par)',/,
     +          '       4 = symmetric (par) probability interval',/,
     +          '       5 = conditional variance')

      close(lun)
      return
      end
