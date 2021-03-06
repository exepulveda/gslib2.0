      real function resc(xmin1,xmax1,xmin2,xmax2,x111)
      real*8 rsc
!
! Simple linear rescaling (get a value in coordinate system "2" given
! a value in "1"):
!
      rsc  = dble((xmax2-xmin2)/(xmax1-xmin1))
      resc = xmin2 + real( dble(x111 - xmin1) * rsc )
      return
      end
