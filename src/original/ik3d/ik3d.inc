C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                                                      %
C Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %
C Junior University.  All rights reserved.                             %
C                                                                      %
C The programs in GSLIB are distributed in the hope that they will be  %
C useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
C responsibility to anyone for the consequences of using them or for   %
C whether they serve any particular purpose or work at all, unless he  %
C says so in writing.  Everyone is granted permission to copy, modify  %
C and redistribute the programs in GSLIB, but only under the condition %
C that this notice and the above copyright notice remain intact.       %
C                                                                      %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c-----------------------------------------------------------------------
c
c Parameters:
c
      parameter(MAXNST=4,EPSLON=1.0e-5,VERSION=2.908)
c
c The data and other input variables:
c
      integer  test
c
c
      common /datcom/ nd,nx,ny,nz,xmn,ymn,zmn,xsiz,ysiz,zsiz,
     +                ndmax,ndmin,radius,sang1,sang2,sang3,sanis1,
     +                sanis2,isrot,ivrot,mik,ktype,noct,idbg,
     +                ldbg,iout,lout,tmin,tmax,ivtype,koption,
     +                idhl,idhlj,test
c
c Search variables and data for kriging:
c
      common /srccom/ na,ljack,ixlj,iylj,izlj,ivrlj,nvarij
c
c Variogram Parameters:
c
      common /vargdt/ ncut               
