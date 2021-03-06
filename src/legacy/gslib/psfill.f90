      subroutine psfill(np,x,y,lwidt,gray)
!-----------------------------------------------------------------------
!
!
! CALLING ARGUMENTS:
!
!  x            X location of the center of the box
!  y            Y location of the center of the box
!  np           number of points
!  lwidt        The width of the line (1.0 = dark, 0.5 = light)
!  gray         the grayness of the fill area
!
! NOTES:
!
!  1. The pxmin,pxmax,.. variables are in the standard 1/72 inch 
!     resolution of the postscript page. If a different scale is 
!     going to be used in the printing set pscl to the scale.
!
!
!
!-----------------------------------------------------------------------
      parameter(EPSLON=0.0001)
      real lwidt,x(1),y(1)
!
! Common Block for Postscript Output Unit and Scaling:
!
      common /psdata/ lpsout,pscl,pxmin,pxmax,pymin,pymax,xmin, &
                      xmax,ymin,ymax
!
! Change the line width:
!
      if(pscl.lt.0.01) pscl = 1.0
      width = lwidt/pscl
      write(lpsout,103) width
!
! Start a new path and loop through the points:
!
      write(lpsout,100)
      do i=1,np
            ix = int(resc(xmin,xmax,pxmin,pxmax,x(i))/pscl)
            iy = int(resc(ymin,ymax,pymin,pymax,y(i))/pscl)
            if(i.eq.1) then
                  write(lpsout,101) ix,iy
            else
                  write(lpsout,102) ix,iy
            endif
      end do         
      if(lwidt.le.EPSLON) then
            write(lpsout,104) gray
      else
            write(lpsout,105) gray
      endif
 100  format('n')
 101  format(i5,1x,i5,1x,'m')
 102  format(i5,1x,i5,1x,'l')
 103  format(f6.3,' setlinewidth')
 104  format('c',/,f4.2,' setgray',/,'fill',/,'0.0 setgray')
 105  format('c gsave ',/,f4.2,' setgray',/,'fill',/,'grestore s', &
               /,'0.00 setgray')
      return
      end
