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
      parameter(MAXNST=4,EPSLON=1.0e-20,UNEST=-99.0,VERSION=2.907)
c
c ACORN parameters:
c
      parameter(KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
c
c Static Array Dimensioning:
c
      integer     inoct(8),ltail,middle,utail,sstrat,mults,nmult,ivtype,
     +            test

      real        ltpar,mpar,utpar
c
c Common Blocks:
c
      common /iaco/   ixv(MAXOP1)

      common /grid3d/ nx,ny,nz,nxy,nxyz,xsiz,ysiz,zsiz,xmn,ymn,zmn

      common /data3d/ ivtype,nd,zmin,zmax,ltail,middle,test,
     +                utail,ltpar,mpar,utpar,tmin,tmax

      common /simula/ nsim,lin,lout,idbg,ldbg

      common /srch3d/ radius,radsqd,sang1,sang2,sang3,sanis1,sanis2,
     +                noct,nclose,ndmax,mik,sstrat,nodmax,ncnode,
     +                infoct,inoct,mults,nmult

      common /cova3d/ ncut,ng,isrot,ivrot,cbb

      common /clooku/ nctx,ncty,nctz,nlooku,ncnodex