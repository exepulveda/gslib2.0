!-----------------------------------------------------------------------
!
!                  Kriging of a 3-D Rectangular Grid
!                  *********************************
!
! This file controls both the static dimensioning and allocates most of
! the storage in common blocks.
!
!   MAXSBX    maximum super block nodes in X direction
!   MAXSBY    maximum super block nodes in Y direction
!   MAXSBZ    maximum super block nodes in Z direction
!   MAXDAT    maximum number of data points
!   MAXDT     maximum number of drift terms
!   MAXSAM    maximum number of data points to use in one kriging system
!   MAXDIS    maximum number of discretization points per block
!   MAXNST    maximum number of nested structures
!
!-----------------------------------------------------------------------
!
! User Adjustable Parameters:
!
      parameter(MAXDT  =   9, MAXNST =   4)
!
! Fixed Parameters:
!
      parameter(UNEST = -999.0, EPSLON = 0.000001, MAXROT = MAXNST + 1, &
                VERSION = 2.907)
!
! Static Array Dimensioning:
!
      integer   idrif(MAXDT),nst(1),it(MAXNST),iktype,ncut,test
      real      bv(9),c0(1),cc(MAXNST),aa(MAXNST),ang1(MAXNST), &
                ang2(MAXNST),ang3(MAXNST),anis1(MAXNST),anis2(MAXNST)
      real*8    rotmat(MAXROT,3,3)
!
! The data and other input variables:
!
      common /datcom/ nd,tmin,tmax,nx,ny,nz,xmn,ymn,zmn,test, &
                      xsiz,ysiz,zsiz,ndmax,ndmin,radius,noct,nxdis, &
                      nydis,nzdis,idrif,itrend,ktype,skmean,koption, &
                      idbg,ldbg,iout,lout,iext,lext,iextve,ljack, &
                      idhlj,ixlj,iylj,izlj,ivrlj,iextvj,nvarij
!
! Variogram Parameters:
!
      common /vargdt/ nst,it,c0,cc,aa,ang1,ang2,ang3,anis1,anis2,rotmat
!
! Search variables and data for kriging:
!
      common /srccom/ sang1,sang2,sang3,sanis1,sanis2,isrot,na,ndb, &
                     bv,unbias,iktype,ncut
