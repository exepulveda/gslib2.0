      subroutine strlen(str,MAXLEN,lostr)
!-----------------------------------------------------------------------
!
!      Determine the length of the string minus trailing blanks
!
!
!
!-----------------------------------------------------------------------
      character str*512
      lostr = MAXLEN
      do i=1,MAXLEN
            j = MAXLEN - i + 1
            if(str(j:j).ne.' ') return
            lostr = lostr - 1
      end do
      return
      end
