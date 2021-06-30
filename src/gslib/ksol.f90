      subroutine ksol(nright,neq,nsb,a,r,s,ising)
!-----------------------------------------------------------------------
!
!                Solution of a System of Linear Equations
!                ****************************************
!
!
!
! INPUT VARIABLES:
!
!   nright,nsb       number of columns in right hand side matrix.
!                      for KB2D: nright=1, nsb=1
!   neq              number of equations
!   a()              upper triangular left hand side matrix (stored 
!                      columnwise)
!   r()              right hand side matrix (stored columnwise)
!                      for kb2d, one column per variable
!
!
!
! OUTPUT VARIABLES:
!
!   s()              solution array, same dimension as  r  above.
!   ising            singularity indicator
!                      0,  no singularity problem
!                     -1,  neq .le. 1
!                      k,  a null pivot appeared at the kth iteration
!
!
!
! PROGRAM NOTES:
!
!   1. Requires the upper triangular left hand side matrix.
!   2. Pivots are on the diagonal.
!   3. Does not search for max. element for pivot.
!   4. Several right hand side matrices possible.
!   5. USE for ok and sk only, NOT for UK.
!
!
!-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8   a(*),r(*),s(*)
!
! If there is only one equation then set ising and return:
!
      if(neq.le.1) then
            ising = -1
            return
      endif
!
! Initialize:
!
      tol   = 0.1e-06
      ising = 0
      nn    = neq*(neq+1)/2
      nm    = nsb*neq
      m1    = neq-1
      kk    = 0
!
! Start triangulation:
!
      do k=1,m1
            kk=kk+k
            ak=a(kk)
            if(abs(ak).lt.tol) then
                  ising=k
                  return
            endif
            km1=k-1
            do iv=1,nright
                  nm1=nm*(iv-1)
                  ii=kk+nn*(iv-1)
                  piv=1./a(ii)
                  lp=0
                  do i=k,m1
                        ll=ii
                        ii=ii+i
                        ap=a(ii)*piv
                        lp=lp+1
                        ij=ii-km1
                        do j=i,m1
                              ij=ij+j
                              ll=ll+j
                              a(ij)=a(ij)-ap*a(ll)
                        end do
                        do llb=k,nm,neq
                              in=llb+lp+nm1
                              ll1=llb+nm1
                              r(in)=r(in)-ap*r(ll1)
                        end do
                  end do
            end do
      end do
!
! Error checking - singular matrix:
!
      ijm=ij-nn*(nright-1)
      if(abs(a(ijm)).lt.tol) then
            ising=neq
            return
      endif
!
! Finished triangulation, start solving back:
!
      do iv=1,nright
            nm1=nm*(iv-1)
            ij=ijm+nn*(iv-1)
            piv=1./a(ij)
            do llb=neq,nm,neq
                  ll1=llb+nm1
                  s(ll1)=r(ll1)*piv
            end do
            i=neq
            kk=ij
            do ii=1,m1
                  kk=kk-i
                  piv=1./a(kk)
                  i=i-1
                  do llb=i,nm,neq
                        ll1=llb+nm1
                        in=ll1
                        ap=r(in)
                        ij=kk
                        do j=i,m1
                              ij=ij+j
                              in=in+1
                              ap=ap-a(ij)*s(in)
                        end do
                        s(ll1)=ap*piv
                  end do
            end do
      end do
!
! Finished solving back, return:
!
      return
      end
