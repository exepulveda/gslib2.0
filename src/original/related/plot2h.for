      program main
c-----------------------------------------------------------------------
c
c  DOS Program to perform what the "plot" script does on a UNIX machine
c
c
c-----------------------------------------------------------------------
      use msflib
      
      parameter(MAXLEN=1024)
      data      lin/2/,lout/3/
      character str*1024,str2*80,str3*80
      integer   nsize
      
      call getarg(1,str)
      call getarg(2,str2)
      call getarg(3,str3)

      open(lin,file=str,status='OLD')
      open(lout,file=str3)
      write(lout,*)'%!'
      write(lout,*)
     +    'gsave 1.0 72 mul 4.5 72 mul translate 0.75 0.75 scale'
      nsize = 0
 80   read(lin,*,end=90,err=99)
      nsize = nsize + 1
      go to 80
 90   continue
      rewind(lin)
      nsize = nsize -2
      do i = 1,11
            read(lin,*)
      end do
      do i = 12, nsize
            read(lin,'(a)',err=99)str
            call strlen(str,MAXLEN,lostr)
            write(lout,'(a)')str(1:lostr)
      end do
      write(lout,*)'grestore'
      close(lin)
      
      open(lin,file=str2,status='OLD')
      write(lout,*)
     +    'gsave 4.5 72 mul 4.5 72 mul translate 0.75 0.75 scale'
      nsize = 0
 81   read(lin,*,end=91,err=99)
      nsize = nsize + 1
      go to 81
 91   continue
      rewind(lin)
      nsize = nsize -2
      do i = 1,11
            read(lin,*)
      end do
      do i = 12, nsize
            read(lin,'(a)',err=99)str
            call strlen(str,MAXLEN,lostr)
            write(lout,'(a)')str(1:lostr)
      end do
      close(lin)
      
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
