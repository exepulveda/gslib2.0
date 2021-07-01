      subroutine chknam(str,len)
!-----------------------------------------------------------------------
!
!                   Check for a Valid File Name
!                   ***************************
!
! This subroutine takes the character string "str" of length "len" and
! removes all leading blanks and blanks out all characters after the
! first blank found in the string (leading blanks are removed first).
!
!
!
!-----------------------------------------------------------------------
      parameter (MAXLEN=512)
      character str(MAXLEN)*1
!
! Find first two blanks and blank out remaining characters:
!
      do i=1,len-1
            if(str(i)  .eq.' '.and. &
               str(i+1).eq.' ') then
                  do j=i+1,len
                        str(j) = ' '
                  end do
                  go to 2
            end if
      end do
 2    continue
!
! Look for "-fi" for file
!
      do i=1,len-2
            if(str(i)  .eq.'-'.and. &
               str(i+1).eq.'f'.and. &
               str(i+2).eq.'i') then
                  do j=i+1,len
                        str(j) = ' '
                  end do
                  go to 3
            end if
      end do
 3    continue
!
! Look for "\fi" for file
!
      do i=1,len-2
            if(str(i)  .eq.'\'.and. &
               str(i+1).eq.'f'.and. &
               str(i+2).eq.'i') then
                  do j=i+1,len
                        str(j) = ' '
                  end do
                  go to 4
            end if
      end do
 4    continue
!
! Return with modified file name:
!
      return
      end
