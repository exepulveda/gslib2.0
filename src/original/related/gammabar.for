      program gammabar
c-----------------------------------------------------------------------
c
c                  Average Variogram for a Volume
c                  ******************************
c
c This program calculates the gammabar value from a semivariogram model
c
c
c
c
c
c-----------------------------------------------------------------------
      parameter(MAXNST=4, MDIR=10, DEG2RAD=3.14159265/180.0, 
     +          MAXROT=MAXNST+1, EPSLON = 1.0e-20, VERSION=2.905)

      real*8     rotmat(MAXROT,3,3)
      real       aa(MAXNST),c0(1),cc(MAXNST),ang1(MAXNST),ang2(MAXNST),
     +           ang3(MAXNST),anis1(MAXNST),anis2(MAXNST),maxcov
      integer    nst(1),it(MAXNST)
      character  str*40
      logical    testfl
      data       lin/1/
c
c Note VERSION number:
c
      write(*,9999) VERSION
 9999 format(/' GAMMABAR Version: ',f5.3/)
c
c Get the name of the parameter file - try the default name if no input:
c
      write(*,*) 'Which parameter file do you want to use?'
      read (*,'(a20)') str(1:20)
      if(str(1:1).eq.' ') str(1:20) = 'gammabar.par        '
      inquire(file=str(1:20),exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'gammabar.par        ') then
                  write(*,*) '        creating a blank parameter file'
                  call makepar
                  write(*,*)
            end if
            stop
      endif
      open(lin,file=str(1:20),status='OLD')
c
c Find Start of Parameters:
c
 1    read(lin,'(a4)',end=98) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
c
c Read Input Parameters:
c

      read(lin,*,err=98) xsiz,ysiz,zsiz
      write(*,*) ' size of block = ',xsiz,ysiz,zsiz

      read(lin,*,err=98) nx,ny,nz
      write(*,*) ' discretization= ',nx,ny,nz

      nx = max(nx,1)
      ny = max(ny,1)
      nz = max(nz,1)

      read(lin,*,err=98) nst(1),c0(1)
      write(*,*) ' nst, c0 = ',nst(1),c0(1)
      if(nst(1).le.0) then
            write(*,9997) nst(1)
 9997       format(' nst must be at least 1, it has been set to ',i4,/,
     +             ' The c or a values can be set to zero')
            stop
      endif
      do i=1,nst(1)
            read(lin,*,err=98) it(i),cc(i),ang1(i),ang2(i),ang3(i)
            read(lin,*,err=98) aa(i),aa1,aa2
            anis1(i) = aa1 / max(aa(i),EPSLON)
            anis2(i) = aa2 / max(aa(i),EPSLON)
            if(it(i).eq.4) then
                  if(aa(i).lt.0.0) stop ' INVALID power variogram'
                  if(aa(i).gt.2.0) stop ' INVALID power variogram'
            end if
            write(*,*) ' it,cc,ang[1,2,3] =  ',it(i),cc(i),
     +                   ang1(i),ang2(i),ang3(i)
            write(*,*) ' a1 a2 a3 = ',aa(i),aa1,aa2
            write(*,*) ' anis 1 2 = ',anis1(i),anis2(i)
      end do
      close(lin)
c
c Get ready to go:
c
      do is=1,nst(1)
            call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is),
     +                  is,MAXROT,rotmat)
      end do
      call cova3(0.0,0.0,0.0,0.0,0.0,0.0,1,nst,MAXNST,c0,it,cc,aa,
     +                 1,MAXROT,rotmat,cmax,maxcov)
c
c Discretization parameters:
c
      xsz   = xsiz / real(nx)
      xmn   = xsz  / 2.0
      xzero = xsz  * 0.0001
      ysz   = ysiz / real(ny)
      ymn   = ysz  / 2.0
      yzero = ysz  * 0.0001
      zsz   = zsiz / real(nz)
      zmn   = zsz  / 2.0
      zzero = zsz  * 0.0001

      write(*,*)
      write(*,*) 'X direction size ',xsiz,' number ',nx,' offset ',xzero
      write(*,*) 'Y direction size ',ysiz,' number ',ny,' offset ',yzero
      write(*,*) 'Z direction size ',zsiz,' number ',nz,' offset ',zzero
c
c Calculate gammabar:
c
      gammab = 0.0
      do ix=1,nx
      xxi = xmn + real(ix-1)*xsz+xzero
      do iy=1,ny
      yyi = ymn + real(iy-1)*ysz+yzero
      do iz=1,nz
      zzi = zmn + real(iz-1)*zsz+zzero

            do jx=1,nx
            xxj = xmn + real(jx-1)*xsz
            do jy=1,ny
            yyj = ymn + real(jy-1)*ysz
            do jz=1,nz
            zzj = zmn + real(jz-1)*zsz

            call cova3(xxi,yyi,zzi,xxj,yyj,zzj,1,nst,MAXNST,c0,it,
     +                 cc,aa,1,MAXROT,rotmat,cmax,cov)
            gammab = gammab + (maxcov - cov)

            end do
            end do
            end do
      end do
      end do
      end do

      gammab = gammab / (real(nx*ny*nz)**2)

      write(*,*)
      write(*,*) ' Average variogram = ',gammab
c
c Finished:
c
      write(*,9998) VERSION
 9998 format(/' GAMMABAR Version: ',f5.3, ' Finished'/)
      stop
 98   stop 'Error in parameter file somewhere'
      end



      subroutine cova3(x1,y1,z1,x2,y2,z2,ivarg,nst,MAXNST,c0,it,cc,aa,
     +                 irot,MAXROT,rotmat,cmax,cova)
c-----------------------------------------------------------------------
c
c                    Covariance Between Two Points
c                    *****************************
c
c This subroutine calculated the covariance associated with a variogram
c model specified by a nugget effect and nested varigoram structures.
c The anisotropy definition can be different for each nested structure.
c
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         coordinates of first point
c   x2,y2,z2         coordinates of second point
c   nst(ivarg)       number of nested structures (maximum of 4)
c   ivarg            variogram number (set to 1 unless doing cokriging
c                       or indicator kriging)
c   MAXNST           size of variogram parameter arrays
c   c0(ivarg)        isotropic nugget constant
c   it(i)            type of each nested structure:
c                      1. spherical model of range a;
c                      2. exponential model of parameter a;
c                           i.e. practical range is 3a
c                      3. gaussian model of parameter a;
c                           i.e. practical range is a*sqrt(3)
c                      4. power model of power a (a must be gt. 0  and
c                           lt. 2).  if linear model, a=1,c=slope.
c   cc(i)            multiplicative factor of each nested structure.
c                      (sill-c0) for spherical, exponential,and gaussian
c                      slope for linear model.
c   aa(i)            parameter "a" of each nested structure.
c   irot             index of the rotation matrix for the first nested 
c                    structure (the second nested structure will use
c                    irot+1, the third irot+2, and so on)
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c
c
c OUTPUT VARIABLES:
c
c   cmax             maximum covariance
c   cova             covariance between (x1,y1,z1) and (x2,y2,z2)
c
c
c
c EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
c                      rotmat    computes rotation matrix for distance
c-----------------------------------------------------------------------
      parameter(PI=3.14159265,PMX=1.e10,EPSLON=1.e-10)
      integer   nst(*),it(*)
      real      c0(*),cc(*),aa(*)
      real*8    rotmat(MAXROT,3,3),hsqd,sqdist
c
c Calculate the maximum covariance value (used for zero distances and
c for power model covariance):
c
      istart = 1 + (ivarg-1)*MAXNST
      cmax   = c0(ivarg)
      do is=1,nst(ivarg)
            ist = istart + is - 1
            if(it(ist).eq.4) then
                  cmax = cmax + PMX
            else
                  cmax = cmax + cc(ist)
            endif
      end do
c
c Check for "zero" distance, return with cmax if so:
c
      hsqd = sqdist(x1,y1,z1,x2,y2,z2,irot,MAXROT,rotmat)
      if(real(hsqd).lt.EPSLON) then
            cova = cmax
            return
      endif
c
c Loop over all the structures:
c
      cova = 0.0
      do is=1,nst(ivarg)
            ist = istart + is - 1
c
c Compute the appropriate distance:
c
            if(ist.ne.1) then
                  ir = min((irot+is-1),MAXROT)
                  hsqd=sqdist(x1,y1,z1,x2,y2,z2,ir,MAXROT,rotmat)
            end if
            h = real(dsqrt(hsqd))
c
c Spherical Variogram Model?
c
            if(it(ist).eq.1) then
                  hr = h/aa(ist)
                  if(hr.lt.1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
c
c Exponential Variogram Model?
c
            else if(it(ist).eq.2) then
                  cova = cova + cc(ist)*exp(-3.0*h/aa(ist))
c
c Gaussian Variogram Model?
c
            else if(it(ist).eq.3) then
                  cova = cova + cc(ist)*exp(-(3.0*h/aa(ist))
     +                                      *(3.0*h/aa(ist)))
c
c Power Variogram Model?
c
            else if(it(ist).eq.4) then
                  cova = cova + cmax - cc(ist)*(h**aa(ist))
c
c Hole Effect Model?
c
            else if(it(ist).eq.5) then
c                 d = 10.0 * aa(ist)
c                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
                  cova = cova + cc(ist)*cos(h/aa(ist)*PI)
            endif
      end do
c
c Finished:
c
      return
      end



      real*8 function sqdist(x1,y1,z1,x2,y2,z2,ind,MAXROT,rotmat)
c-----------------------------------------------------------------------
c
c    Squared Anisotropic Distance Calculation Given Matrix Indicator
c    ***************************************************************
c
c This routine calculates the anisotropic distance between two points
c  given the coordinates of each point and a definition of the
c  anisotropy.
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         Coordinates of first point
c   x2,y2,z2         Coordinates of second point
c   ind              The rotation matrix to use
c   MAXROT           The maximum number of rotation matrices dimensioned
c   rotmat           The rotation matrices
c
c
c
c OUTPUT VARIABLES:
c
c   sqdist           The squared distance accounting for the anisotropy
c                      and the rotation of coordinates (if any).
c
c
c NO EXTERNAL REFERENCES
c
c
c Author: C. Deutsch                                Date: September 1989
c-----------------------------------------------------------------------
      real*8 rotmat(MAXROT,3,3),cont,dx,dy,dz
c
c Compute component distance vectors and the squared distance:
c
      dx = dble(x1 - x2)
      dy = dble(y1 - y2)
      dz = dble(z1 - z2)
      sqdist = 0.0
      do i=1,3
            cont   = rotmat(ind,i,1) * dx
     +             + rotmat(ind,i,2) * dy
     +             + rotmat(ind,i,3) * dz
            sqdist = sqdist + cont * cont
      end do
      return
      end



      subroutine setrot(ang1,ang2,ang3,anis1,anis2,ind,MAXROT,rotmat)
c-----------------------------------------------------------------------
c
c              Sets up an Anisotropic Rotation Matrix
c              **************************************
c
c Sets up the matrix to transform cartesian coordinates to coordinates
c accounting for angles and anisotropy (see manual for a detailed
c definition):
c
c
c INPUT PARAMETERS:
c
c   ang1             Azimuth angle for principal direction
c   ang2             Dip angle for principal direction
c   ang3             Third rotation angle
c   anis1            First anisotropy ratio
c   anis2            Second anisotropy ratio
c   ind              matrix indicator to initialize
c   MAXROT           maximum number of rotation matrices dimensioned
c   rotmat           rotation matrices
c
c
c NO EXTERNAL REFERENCES
c
c
c Author: C. Deutsch                                Date: September 1989
c-----------------------------------------------------------------------
      parameter(DEG2RAD=3.141592654/180.0,EPSLON=1.e-20)
      real*8    rotmat(MAXROT,3,3),afac1,afac2,sina,sinb,sint,
     +          cosa,cosb,cost
c
c Converts the input angles to three angles which make more
c  mathematical sense:
c
c         alpha   angle between the major axis of anisotropy and the
c                 E-W axis. Note: Counter clockwise is positive.
c         beta    angle between major axis and the horizontal plane.
c                 (The dip of the ellipsoid measured positive down)
c         theta   Angle of rotation of minor axis about the major axis
c                 of the ellipsoid.
c
      if(ang1.ge.0.0.and.ang1.lt.270.0) then
            alpha = (90.0   - ang1) * DEG2RAD
      else
            alpha = (450.0  - ang1) * DEG2RAD
      endif
      beta  = -1.0 * ang2 * DEG2RAD
      theta =        ang3 * DEG2RAD
c
c Get the required sines and cosines:
c
      sina  = dble(sin(alpha))
      sinb  = dble(sin(beta))
      sint  = dble(sin(theta))
      cosa  = dble(cos(alpha))
      cosb  = dble(cos(beta))
      cost  = dble(cos(theta))
c
c Construct the rotation matrix in the required memory:
c
      afac1 = 1.0 / dble(max(anis1,EPSLON))
      afac2 = 1.0 / dble(max(anis2,EPSLON))
      rotmat(ind,1,1) =       (cosb * cosa)
      rotmat(ind,1,2) =       (cosb * sina)
      rotmat(ind,1,3) =       (-sinb)
      rotmat(ind,2,1) = afac1*(-cost*sina + sint*sinb*cosa)
      rotmat(ind,2,2) = afac1*(cost*cosa + sint*sinb*sina)
      rotmat(ind,2,3) = afac1*( sint * cosb)
      rotmat(ind,3,1) = afac2*(sint*sina + cost*sinb*cosa)
      rotmat(ind,3,2) = afac2*(-sint*cosa + cost*sinb*sina)
      rotmat(ind,3,3) = afac2*(cost * cosb)
c
c Return to calling program:
c
      return
      end



      subroutine makepar
c-----------------------------------------------------------------------
c
c                      Write a Parameter File
c                      **********************
c
c
c
c-----------------------------------------------------------------------
      lun = 99
      open(lun,file='gammabar.par')
      write(lun,10)
 10   format('                  Parameters for GAMMABAR',/,
     +       '                  ***********************',/,/,
     +       'START OF PARAMETERS:')

      write(lun,11)
 11   format('160.0  160.0   50.0           ',
     +       '-X,Y,Z size of block')
      write(lun,12)
 12   format(' 21     21     21             ',
     +       '-X,Y,Z discretization')
      write(lun,13)
 13   format('2   0.0                       ',
     +       '-nst, nugget effect')
      write(lun,14)
 14   format('1 0.6    0.0     0.0   0.0    ',
     +       '-  it,cc,ang1,ang2,ang3')
      write(lun,15)
 15   format('      1000.0  3000.0  12.0    ',
     +       '-  a_hmax, a_hmin, a_vert')
      write(lun,16)
 16   format('1 0.4    0.0     0.0   0.0    ',
     +       '-  it,cc,ang1,ang2,ang3')
      write(lun,17)
 17   format('     30000.0  6000.0  50.0    ',
     +       '-  a_hmax, a_hmin, a_vert')

      close(lun)
      return
      end
