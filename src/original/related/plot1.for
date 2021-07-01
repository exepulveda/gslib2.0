      program main
c-----------------------------------------------------------------------
c
c  DOS Program to perform what the "plot" script does on a UNIX machine
c
c
c-----------------------------------------------------------------------
      use msflib
      
      parameter(MAXLEN=1024)
      data      lin/1/,lout/2/
      character str*1024,str2*80
      integer   nsize
      
      call getarg(1,str)
      call getarg(2,str2)
c
c Header on output file:
c
      open(lout,file=str2)
      write(lout,100)'%!'
      write(lout,*)
 100  format('%!',/,'gsave 1.25 72 mul 3.0 72 mul translate',/,
     +       '1.5 1.5 scale')
c
c Keep only what is required:
c
      open(lin,file=str,status='OLD')
      nsize = 0
 8    read(lin,*,end=9,err=99)
      nsize = nsize + 1
      go to 8
 9    continue
      rewind(lin)

      do i = 1,11
            read(lin,*)
      end do
      do i = 12, nsize-2
            read(lin,'(a)',err=99)str
            call strlen(str,MAXLEN,lostr)
            write(lout,'(a)')str(1:lostr)
      end do

      write(lout,101)
 101  format('grestore',/,'%%%EOF',/,'showpage')

      stop
 99   write(*,*)'ERROR in Postscript File.'
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
