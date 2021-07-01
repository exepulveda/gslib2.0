      subroutine dlocate(xx,n,is,ie,x,j)
!-----------------------------------------------------------------------
!
! Given an array "xx" of length "n", and given a value "x", this routine
! returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
! must be monotonic, either increasing or decreasing.  j=0 or j=n is
! returned to indicate that x is out of range.
!
! Modified to set the start and end points by "is" and "ie" 
!
! Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90.
!-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension xx(n)
!
! Initialize lower and upper methods:
!
      jl = is-1
      ju = ie
!
! If we are not done then compute a midpoint:
!
 10   if(ju-jl.gt.1) then
            jm = (ju+jl)/2
!
! Replace the lower or upper limit with the midpoint:
!
            if((xx(ie).gt.xx(is)).eqv.(x.gt.xx(jm))) then
                  jl = jm
            else
                  ju = jm
            endif
            go to 10
      endif
!
! Return with the array index:
!
      j = jl
      return
      end
