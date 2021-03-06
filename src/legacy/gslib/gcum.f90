      real function gcum(x)
!-----------------------------------------------------------------------
!
! Evaluate the standard normal cdf given a normal deviate x.  gcum is
! the area under a unit normal curve to the left of x.  The results are
! accurate only to about 5 decimal places.
!
!
!-----------------------------------------------------------------------
      z = x
      if(z.lt.0.) z = -z
      t    = 1./(1.+ 0.2316419*z)
      gcum = t*(0.31938153   + t*(-0.356563782 + t*(1.781477937 + &
             t*(-1.821255978 + t*1.330274429))))
      e2   = 0.
!
!  6 standard deviations out gets treated as infinity:
!
      if(z.le.6.) e2 = exp(-z*z/2.)*0.3989422803
      gcum = 1.0- e2 * gcum
      if(x.ge.0.) return
      gcum = 1.0 - gcum
      return
      end
