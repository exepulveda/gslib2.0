      program main
c-----------------------------------------------------------------------
c
c         GETDATA: Extract Well Data Within Stratigraphic Layer
c         *****************************************************
c
c
c
c
c
c
c AUTHOR: Clayton V. Deutsch                                  DATE: 1999
c-----------------------------------------------------------------------
      use       msflib
      parameter(MAXLEN=132,MAXWELL=1000,VERSION=1.000)
c
c Dimensioning:
c
      real      var(200),zet(MAXWELL),zeb(MAXWELL),zrt(MAXWELL),
     +          zrb(MAXWELL)
      integer   zwid(MAXWELL)
      character datafl*1024,outfl*1024,str*1024
      logical   testfl
      lin  = 1
      lout = 2
c
c "blank-out" file name character strings:
c
      do i=1,1024
            datafl(i:i) =' '
            outfl(i:i)  =' '
            str(i:i)    =' '
      end do
c
c Note VERSION number before anything else:
c
      write(*,9999) VERSION
 9999 format(/' GETDATA Version: ',f5.3/)
c
c Get the name of the parameter file - try the default name if no input:
c
      call getarg(1,str)
      if(str(1:1).eq.' ')then
            write(*,*) 'Which parameter file do you want to use?'
            read (*,'(a)') str
      end if
      if(str(1:1).eq.' ')str='getdata.par                             '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'getdata.par         ') then
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
 1    read(lin,'(a4)',end=98) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
c
c Read Input Parameters:
c
      read(lin,'(a)',err=98) datafl
      call chknam(datafl,1024)
      write(*,*) ' data file = ',datafl(1:40)

      read(lin,*,err=98) iwell,idepth
      write(*,*) ' columns = ',iwell,idepth

      read(lin,*,err=98) tmin,tmax
      write(*,*) ' trimming limits = ',tmin,tmax

      read(lin,'(a)',err=98) outfl
      call chknam(outfl,1024)
      write(*,*) ' file for output = ',outfl(1:40)

      read(lin,*,err=98) nwell
      write(*,*) ' number of wells = ',nwell
      if(nwell.gt.MAXWELL) stop ' too many wells!'

      do i=1,nwell
            read(lin,*,err=98) zwid(i),zet(i),zeb(i),zrt(i),zrb(i)
            write(*,*) ' well = ',zwid(i)
      end do

      close(lin)
c
c Make sure the file exists:
c
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR: ',datafl(1:40),' does not exist'
            write(*,*) '       this file is needed! '
            stop
      endif
c
c Open up the file with input data:
c
      open(lin,file=datafl,status='UNKNOWN')
      open(lout,file=outfl,status='UNKNOWN')
      read(lin,'(a40)',err=99) str(1:40)
      write(lout,100)          str(1:40)
 100  format('GETDATA:',a40)
      read(lin,*,err=99)       nvari
      write(lout,'(i3)')       nvari+1
      do i=1,nvari
            read(lin,'(a40)',err=99) str(1:40)
            write(lout,'(a40)')      str(1:40)
      end do
      write(lout,101)
 101  format('Layer-Specific Depth')
c
c Read and write all the data until the end of the file:
c
      nin  = 0
      nout = 0
 7    read(lin,*,end=9) (var(i),i=1,nvari)
c
c Check well ID and see if in or out:
c
      if(var(idepth).lt.tmin.or.var(idepth).ge.tmax) then
            nout = nout + 1
            go to 7
      endif
c
c Find out if this data belongs to a recognized well:
c
      depth = var(idepth)
      kwell = int(var(iwell)+0.5)
      do i=1,nwell
            if(zwid(i).eq.kwell) then
                  if((depth.lt.zet(i).and.depth.lt.zeb(i)).or.
     +               (depth.gt.zet(i).and.depth.gt.zeb(i))) then
                        nout = nout + 1
                        go to 7	
                  end if
                  zrel = abs((depth-zrb(i))/(zrt(i)-zrb(i)))
                  backspace lin
                  read(lin,'(a)') str
                  call strlen(str,MAXLEN,lostr)
                  write(lout,'(a,1x,f12.5)') str(1:lostr),zrel
                  nin = nin + 1
                  go to 7
            end if
      end do
      nout = nout + 1
      go to 7
 9    continue
      write(*,*) ' There were ',nin,' data in the layer'
      write(*,*) ' There were ',nout,' data outside the layer'
c
c Finished:
c
      close(lin)
      close(lout)
      write(*,9998) VERSION
 9998 format(/' GETDATA Version: ',f5.3, ' Finished'/)
      stop
 98   stop ' ERROR in parameter file'
 99   stop ' ERROR in data file'
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
      open(lun,file='nscore.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for GETDATA',/,
     +       '                  **********************',/,/,
     +       'START OF PARAMETERS:')

      write(lun,11)
 11   format('sample.dat               ',
     +       '-file with data')
      write(lun,12)
 12   format('1   2                    ',
     +       '-  columns for well ID and depth')
      write(lun,13)
 13   format('-1.0e21   1.0e21         ',
     +       '-  trimming limits')
      write(lun,14)
 14   format('getdata.out              ',
     +       'file for layer data')
      write(lun,16)
 16   format('4                        ',
     +       '-number of well IDs to process')
      write(lun,17)
 17   format('100 0.0 100.0 0.0 100.0  ',
     +       '-  well, existing top/base, restored top/base')
      write(lun,18)
 18   format('101 0.0 100.0 0.0 100.0  ',
     +       '-  well, existing top/base, restored top/base')
      write(lun,19)
 19   format('102 0.0 100.0 0.0 100.0  ',
     +       '-  well, existing top/base, restored top/base')
      write(lun,20)
 20   format('103 0.0 100.0 0.0 100.0  ',
     +       '-  well, existing top/base, restored top/base')

      close(lun)
      return
      end



      subroutine strlen(str,MAXLEN,lostr)
c-----------------------------------------------------------------------
c
c      Determine the length of the string minus trailing blanks
c
c
c
c-----------------------------------------------------------------------
      character str*1024
      lostr = MAXLEN
      do i=1,MAXLEN
            j = MAXLEN - i + 1
            if(str(j:j).ne.' ') return
            lostr = lostr - 1
      end do
      return
      end



      subroutine chknam(str,len)
c-----------------------------------------------------------------------
c
c                   Check for a Valid File Name
c                   ***************************
c
c This subroutine takes the character string "str" of length "len" and
c attempts to blank out all characters that relate to a comment in a
c GSLIB parameter file.  The modifications were made to accomodate
c long file names with blanks that happen in Windows.
c
c
c
c-----------------------------------------------------------------------
      character str*1024
c
c Find "-fil" and blank out the remaining characters:
c
      do i=1,len-3
            if(str(i:i+3).eq.'-fil') then
                  do j=i,len
                        str(j:j) = ' '
                  end do
                  go to 1
            end if
      end do
 1    continue
c
c Find first two blanks and blank out the remaining characters:
c
      do i=1,len-1
            if(str(i:i+1).eq.'  ') then
                  do j=i,len
                        str(j:j) = ' '
                  end do
                  go to 2
            end if
      end do
 2    continue
c
c Return with modified file name:
c
      return
      end
