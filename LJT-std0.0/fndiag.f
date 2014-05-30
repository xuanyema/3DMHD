 	subroutine fndiag
c       This subroutine calculates the `force normal' over the 
c       system.  It is based upon the subroutine diag1.f of the
c       3D MHD code.
c       Version 15p
c       This version is based on Version 14p.  Version 15p differs from
c       its predecessor in the following ways:
c       1) The value of ixx has been changed.
c       2) The dimensionality of fnintg is different; it does not
c          include the ldiag dimension.  Furthermore, this variable
c          is defined in this subroutine, rather than in the misflin
c          header file.
c       3) The limits of the loops used in calculating the integrand
c          and integral of the force norm have been changed.
c       17 June 2003
c***********************************************************************
	include 'misflin'
c
        integer  iz,ix,iy,nxe,icomp,ixx
	integer  nyhalb,trgtx,trgty
        real     bmax(4,4), bmin(4,4), vmax(4,4), vmin(4,4),
     +           jmax(4,4), jmin(4,4), emax(4,4), emin(4,4),
     +           jemin(2,4),jemax(2,4),prho(4,4),
     +           fmax(nx),fmin(nx)
	real     tmhdfx(nx,ny,nz),tmhdfy(nx,ny,nz),
     +           tmhdfz(nx,ny,nz)
	real     fnintg(nx,ny,nz)
	real     slmax

c.....................................................................
c   Compute the differential volume elements (or, rather, their 
c   practical analogues here) over the grid.  Take special care at 
c   the boundaries.
c.....................................................................
c     ixx = (nx+1)/2
      ixx = 3
      nxe = (nx+1)/2
      nyhalb = (ny+1)/2
      trgtx = nx - 1
      trgty = nyhalb
c     write(6,*) 'Target diagnostic line defined by'
c     write(6,247) trgtx,trgty
c     write(6,248) x(trgtx),y(trgty)
      do 20 iz = 1,nz
      do 20 iy = 1,ny
      do 20 ix = 1,nx
	hilf(ix,iy,iz) = 0.125/difx(ix)/dify(iy)/difz(iz)
	if (ix .eq. 2 .or. ix .eq. nx1) 
     +          hilf(ix,iy,iz) = 0.5*hilf(ix,iy,iz) 
	if (iy .eq. 2 .or. iy .eq. ny1) 
     +          hilf(ix,iy,iz) = 0.5*hilf(ix,iy,iz) 
	if (iz .eq. 2 .or. iz .eq. nz1) 
     +          hilf(ix,iy,iz) = 0.5*hilf(ix,iy,iz) 
  20  continue   

c......................................................................
c   Compute the components of the current density.
c......................................................................
      do 25 iz = 1,nz
      do 25 iy = 1,ny
      do 25 ix = 1,nx	  
         vx(ix,iy,iz) = sx(ix,iy,iz)*hilf(ix,iy,iz)
         vy(ix,iy,iz) = sy(ix,iy,iz)*hilf(ix,iy,iz)
         vz(ix,iy,iz) = sz(ix,iy,iz)*hilf(ix,iy,iz)
  25  continue
      do 30 iz = 1,nzz
      do 30 iy = 1,ny
      do 30 ix = 1,nx	  
         help(ix,iy,iz) = u(ix,iy,iz)**gamma + 0.5*(
     +                        + bx(ix,iy,iz)*bx(ix,iy,iz)
     +                        + by(ix,iy,iz)*by(ix,iy,iz)
     +                        + bz(ix,iy,iz)*bz(ix,iy,iz) ) 
  30  continue
      do 35 iz = 1,nz
      do 35 iy = 1,ny
      do 35 ix = 2,nx
	 hilfx(ix,iy,iz) = help(ix,iy,iz)
     +                   - bx(ix,iy,iz)*bx(ix,iy,iz)
	 hilfy(ix,iy,iz) = help(ix,iy,iz)
     +                   - by(ix,iy,iz)*by(ix,iy,iz)
	 hilfz(ix,iy,iz) = help(ix,iy,iz)
     +                   - bz(ix,iy,iz)*bz(ix,iy,iz)
	 feldx(ix,iy,iz) = - by(ix,iy,iz)*bz(ix,iy,iz)
	 feldy(ix,iy,iz) = - bz(ix,iy,iz)*bx(ix,iy,iz)
	 feldz(ix,iy,iz) = - bx(ix,iy,iz)*by(ix,iy,iz)
  35  continue
      do 40 iz = 2,nz1
      do 40 iy = 2,ny1
      do 40 ix = 2,nx1
	 tmhdfx(ix,iy,iz)=difx(ix)*(hilfx(ix+1,iy,iz)-hilfx(ix-1,iy,iz))
     +                 + dify(iy)*(feldz(ix,iy+1,iz)-feldz(ix,iy-1,iz))
     +                 + difz(iz)*(feldy(ix,iy,iz+1)-feldy(ix,iy,iz-1))
	 tmhdfy(ix,iy,iz)=dify(iy)*(hilfy(ix,iy+1,iz)-hilfy(ix,iy-1,iz))
     +                 + difz(iz)*(feldx(ix,iy,iz+1)-feldx(ix,iy,iz-1))
     +                 + difx(ix)*(feldz(ix+1,iy,iz)-feldz(ix-1,iy,iz))
         tmhdfz(ix,iy,iz)=difz(iz)*(hilfz(ix,iy,iz+1)-hilfz(ix,iy,iz-1))
     +                 + difx(ix)*(feldy(ix+1,iy,iz)-feldy(ix-1,iy,iz))
     +                 + dify(iy)*(feldx(ix,iy+1,iz)-feldx(ix,iy-1,iz))
  40  continue


      do 50 iz = 2,nz1
      do 50 iy = 2,ny1
      do 50 ix = 2,nx1
        feldx(ix,iy,iz) = dify(iy)*(bz(ix,iy+1,iz)-bz(ix,iy-1,iz))
     +                     -difz(iz)*(by(ix,iy,iz+1)-by(ix,iy,iz-1))
        feldy(ix,iy,iz) = difz(iz)*(bx(ix,iy,iz+1)-bx(ix,iy,iz-1))
     +                     -difx(ix)*(bz(ix+1,iy,iz)-bz(ix-1,iy,iz))
        feldz(ix,iy,iz) = difx(ix)*(by(ix+1,iy,iz)-by(ix-1,iy,iz))
     +                     -dify(iy)*(bx(ix,iy+1,iz)-bx(ix,iy-1,iz))
c
  50  continue
      do 55 iz = 1,nz
      do 55 iy = 1,ny
      do 55 ix = 1,nx	  
         help(ix,iy,iz) = u(ix,iy,iz)**gamma 
  55  continue

c......................................................................
c   Note concern about factor of 2 for pressure.
c......................................................................

c
c.....................................................................
c   Compute the integrand of the force normal throughout the grid.
c.....................................................................
      do 100 iz = 2,nz-2
      do 100 iy = 3,ny-2
      do 100 ix = ixx,nx-2
	fnintg(ix,iy,iz) =   tmhdfx(ix,iy,iz)*tmhdfx(ix,iy,iz)
     +                     + tmhdfy(ix,iy,iz)*tmhdfy(ix,iy,iz)
     +                     + tmhdfz(ix,iy,iz)*tmhdfz(ix,iy,iz)
  100 continue
c
c
c.....................................................................
c   Integrate over the grid.
c.....................................................................
      do 200 iz = 2,nz-2
      do 200 iy = 3,ny-2
      do 200 ix = ixx,nx-2
	 fnorm(ldiag) = fnorm(ldiag) 
     +                  + hilf(ix,iy,iz)*fnintg(ix,iy,iz)
  200 continue
c
c
c
c......................................................................
c     Determine the maximum of the force norm, as well as the location
c     at which that maximum is attained.
c     Note that the code is based upon that in the subroutine diag1.
c......................................................................
c     Initialize relevant variables.
      fnmax(ldiag) = -100.0
c
c     Identify the maximum value of the force norm density and the
c     location (in terms of the indices) at which that maximum is
c     attained.
      do 220 iz = 2,nz-2
      do 220 iy = 3,ny-2
      do 220 ix = ixx,nx-2
	 if (fnintg(ix,iy,iz).gt.fnmax(ldiag)) then
	    fnmax(ldiag) = fnintg(ix,iy,iz)
	    ixmm(ldiag) = ix
	    iymm(ldiag) = iy
	    izmm(ldiag) = iz
	 endif
 220  continue
c
c     Print the maximum value of the force norm to the standard
c     output, along with the location of that maximum.
      write(6,270) 
     +       fnmax(ldiag),rho(ixmm(ldiag),iymm(ldiag),izmm(ldiag)),time
      write(6,275) ixmm(ldiag),iymm(ldiag),izmm(ldiag)

c

  247 format(1x,'ix = ',i2,2x,'iy = ',i2)
  248 format(1x,'x = ',f14.5,2x,'y = ',f14.5)
  251 format(1x,'ix = ',i2,2x,'iy = ',i2,2x,'iz = ',i2)
  252 format(1x,'x = ',f14.5,2x,'y = ',f14.5,2x,'z = ',f14.5)
  253 format(1x,'difx = ',g14.5,2x,'dify = ',g14.5,2x,'difz = ',g14.5)
  254 format(1x,'bx = ',g14.5,2x,'by = ',g14.5,2x,'bz = ',g14.5)
  255 format(1x,'jx = ',g14.5,2x,'jy = ',g14.5,2x,'jz = ',g14.5)
  256 format(1x,'(F_m)_x = ',g14.5,2x,'(F_m)_y = ',g14.5,
     +       2x,'(F_m)_z= ',g14.5)
  257 format(1x,'(F_p)_x = ',g14.5,2x,'(F_p)_y = ',g14.5,
     +       2x,'(F_p)_z= ',g14.5)
  258 format(1x,'(F_tot)_x = ',g14.5,2x,'(F_tot)_y = ',g14.5,
     +       2x,'(F_tot)_z= ',g14.5)
  259 format((/))
  270 format('Max(fnorm) = ',g14.5,'  Density:',g14.5,'  time:',f10.4)
  275 format(' at (',i5,',',i5,',',i5,')')



      return
      end
