	subroutine grid
c***********************************************************************
	include 'misflin'
c
	integer  ix,ixs,iy,iys,iz,ix0,iy0,iz0
	real     kx,ky,kz,xanf,yanf,zanf,
     +           beta,a1,a2,a3,a4,w,w1,w2
c .................................................................
c     grid computes gridparameters
c         ( x perp to current sheet       -  n koordinate
c           y                             - -m koordinate
c           z anti parallel field         -  l koordinate )
c .................................................................
c
c   x-coordinate
c  (y und z similar) using:
c         x  =  w + a1*sin(beta*w) + a2*sin(2*beta*w)
c         w  =  kx*(ix-ix0)
c         kx =  2.*xmax/(nx-3) ; ix0 = (nx+1)/2
c         beta = pi/xmax
c         a1,a2 determined by grid distance at x=0 (eps(1))
c         difx - coeff 1. derivative
c         ddifx, ddifpx und ddifmx - 2. derivative
c ..............................................................
      write(*,991)
 991  format('vor x gitter')

      kx =  (xmax-xmin)/(nx-3)
      xanf = xmin
      ix0 = 2
      beta = pi/(xmax-xmin)
      if (zentr(1)) then
	xanf = (xmin+xmax)/2.0
	ix0 = (nx+1)/2
	beta = 2*pi/(xmax-xmin)
      end if
      a1 = 0.0
      a2 = 0.0
      a3 = 0.0
      a4 = 0.0
      if (.not. aequi(1)) then
	a1 = 24./17.*(eps(1)-kx)/kx/beta
	a2 =  6./17.*(kx-eps(1))/kx/beta
	a3 =  8./51.*(eps(1)-kx)/kx/beta
	a4 =  3./68.*(kx-eps(1))/kx/beta
      end if
c.......................................
      do 200 ix = 1,nx
         w      =  kx*(ix-ix0)
	 x(ix)  =  w + a1*sin(beta*w) + a2*sin(2*beta*w)
     +               + a3*sin(3.*beta*w) + a4*sin(4.*beta*w)
	 difx(ix) = 0.5/kx/(1. + a1*beta*cos(beta*w)
     +                    + 2.*a2*beta*cos(2.*beta*w)
     +                     + 3.*a3*beta*cos(3.*beta*w)
     +                      + 4.*a4*beta*cos(4.*beta*w) )
         w1 = w + 0.5*kx
	 w2 = w - 0.5*kx
	 ddifpx(ix) = 2.*difx(ix)/kx/(1. + a1*beta*cos(beta*w1)
     +                              + 2.*a2*beta*cos(2.*beta*w1)
     +                               + 3.*a3*beta*cos(3.*beta*w1)
     +                                + 4.*a4*beta*cos(4.*beta*w1) )
	 ddifmx(ix) = 2.*difx(ix)/kx/(1. + a1*beta*cos(beta*w2)
     +                              + 2.*a2*beta*cos(2.*beta*w2)
     +                               + 3.*a3*beta*cos(3.*beta*w2)
     +                                + 4.*a4*beta*cos(4.*beta*w2) )
         ddifx(ix) = .5*(ddifpx(ix)+ddifmx(ix))
  200 continue
c
c  coeff for averages:
c...........................................
      do 250 ix = 2,nx1
	 meanmx(ix) = (x(ix+1)-x(ix))/(x(ix+1)-x(ix-1))
	 meanpx(ix) = (x(ix)-x(ix-1))/(x(ix+1)-x(ix-1))
  250 continue
      meanmx(1)   = meanpx(3)
      meanpx(1)   = meanmx(3)
      meanmx(nx)  = meanpx(nx-2)
      meanpx(nx)  = meanmx(nx-2)
c
c symmetriy for centered nonuniform grid
c.......................................................
      if (zentr(1))   then
	do 260 ix = 1, ix0
	  x(ix)      = ( x(ix)-x(nx+1-ix) )/2.0
	  difx(ix)   = ( difx(ix)+difx(nx+1-ix) )/2.0
	  ddifpx(ix) = ( ddifpx(ix)+ddifmx(nx+1-ix) )/2.0
	  ddifmx(ix) = ( ddifmx(ix)+ddifpx(nx+1-ix) )/2.0
	  ddifx(ix)  = ( ddifx(ix)+ddifx(nx+1-ix) )/2.0
	  meanpx(ix) = ( meanpx(ix)+meanmx(nx+1-ix) )/2.0
	  meanmx(ix) = ( meanmx(ix)+meanpx(nx+1-ix) )/2.0
  260   continue
	  x(ix0) = 0.0
	  ddifpx(ix0) = ( ddifpx(ix0)+ddifmx(ix0) )/2.0
	  ddifmx(ix0) = ( ddifpx(ix0)+ddifmx(ix0) )/2.0
	  ddifx(ix0) =  ( ddifpx(ix0)+ddifmx(ix0) )/2.0
	  meanpx(ix0) = ( meanpx(ix0)+meanmx(ix0) )/2.0
	  meanmx(ix0) = ( meanpx(ix0)+meanmx(ix0) )/2.0
	do 270 ix = ix0+1, nx
	  x(ix)      = -x(nx+1-ix)
	  difx(ix)   = difx(nx+1-ix)
	  ddifpx(ix) = ddifmx(nx+1-ix)
	  ddifmx(ix) = ddifpx(nx+1-ix)
	  ddifx(ix)  = ddifx(nx+1-ix)
	  meanpx(ix) = meanmx(nx+1-ix)
	  meanmx(ix) = meanpx(nx+1-ix)
  270   continue

      end if
      do 300 ix = 1, nx
	x(ix) = xanf + x(ix)
  300 continue
c
c  gridboundary
c..............................................
      x(1)        = 2.*xmin - x(3)
      x(nx)       = 2.*xmax -x(nx-2)
      x(2)        = xmin
      x(nx1)     = xmax
      difx(1)     = difx(3)
      difx(nx)    = difx(nx-2)
      ddifmx(2)   = ( ddifpx(2) + ddifmx(2) ) / 2.0
      ddifmx(nx1)= ( ddifpx(nx1) + ddifmx(nx1) ) / 2.0
      ddifpx(2)   = ddifmx(2)
      ddifpx(nx1)= ddifmx(nx1)
      ddifmx(1)   = ddifpx(3)
      ddifpx(1)   = ddifmx(3)
      ddifmx(nx)  = ddifpx(nx-2)
      ddifpx(nx)  = ddifmx(nx-2)
      meanmx(1)   = meanpx(3)
      meanpx(1)   = meanmx(3)
      meanmx(nx)  = meanpx(nx-2)
      meanpx(nx)  = meanmx(nx-2)
c
      if (zentr(1))   then
	do 360 ix = 1,nx
	  ixs = 1+nx-ix
	  if ( meanpx(ix) .ne. meanmx(ixs) )
     +      write(26,21) ix,meanpx(ix),meanmx(ixs),time
	  if ( ddifpx(ix) .ne. ddifmx(ixs) )
     +      write(26,22) ix,ddifpx(ix),ddifmx(ixs),time
  360   continue
      end if
   21 format(' meanpx(',i3,') =',f12.7,
     +       '   and',f12.7,'   time =',f9.4)
   22 format(' ddifpx(',i3,') =',f12.7,
     +       '   and',f12.7,'   time =',f9.4)

c.......................................................
c   y-coordinate
c      analogue x
c.......................................................
      ky =  (ymax-ymin)/(ny-3)
      yanf = ymin
      iy0 = 2
      beta = pi/(ymax-ymin)
      if (zentr(2)) then
	yanf = (ymin+ymax)/2.0
	iy0 = (ny+1)/2
	beta = 2*pi/(ymax-ymin)
      end if
      a1 = 0.0
      a2 = 0.0
      a3 = 0.0
      a4 = 0.0
      if (.not. aequi(2)) then
	a1 = 24./17.*(eps(2)-ky)/ky/beta
	a2 =  6./17.*(ky-eps(2))/ky/beta
	a3 =  8./51.*(eps(2)-ky)/ky/beta
	a4 =  3./68.*(ky-eps(2))/ky/beta
      end if
c.......................................
      do 400 iy = 1,ny
         w      =  ky*(iy-iy0)
	 y(iy)  =  w + a1*sin(beta*w) + a2*sin(2*beta*w)
     +               + a3*sin(3.*beta*w) + a4*sin(4.*beta*w)
	 dify(iy) = 0.5/ky/(1. + a1*beta*cos(beta*w)
     +                    + 2.*a2*beta*cos(2.*beta*w)
     +                     + 3.*a3*beta*cos(3.*beta*w)
     +                      + 4.*a4*beta*cos(4.*beta*w) )
         w1 = w + 0.5*ky
         w2 = w - 0.5*ky
	 ddifpy(iy) = 2.*dify(iy)/ky/(1. + a1*beta*cos(beta*w1)
     +                              + 2.*a2*beta*cos(2.*beta*w1)
     +                               + 3.*a3*beta*cos(3.*beta*w1)
     +                                + 4.*a4*beta*cos(4.*beta*w1) )
	 ddifmy(iy) = 2.*dify(iy)/ky/(1. + a1*beta*cos(beta*w2)
     +                              + 2.*a2*beta*cos(2.*beta*w2)
     +                               + 3.*a3*beta*cos(3.*beta*w2)
     +                                + 4.*a4*beta*cos(4.*beta*w2) )
	 ddify(iy) = .5*(ddifpy(iy)+ddifmy(iy))
  400 continue
c
c  coeff for averages
c............................................
c
      do 450 iy = 2,ny1
	 meanmy(iy) = (y(iy+1)-y(iy))/(y(iy+1)-y(iy-1))
	 meanpy(iy) = (y(iy)-y(iy-1))/(y(iy+1)-y(iy-1))
  450 continue
      meanmy(1) = meanpy(3)
      meanpy(1) = meanmy(3)
      meanmy(ny)   = meanpy(ny-2)
      meanpy(ny)   = meanmy(ny-2)
c
c centered nonuniform grid
c.......................................................
      if (zentr(2))   then
	do 460 iy = 1, iy0
	  y(iy)      = ( y(iy)-y(ny+1-iy) )/2.0
	  dify(iy)   = ( dify(iy)+dify(ny+1-iy) )/2.0
	  ddifpy(iy) = ( ddifpy(iy)+ddifmy(ny+1-iy) )/2.0
	  ddifmy(iy) = ( ddifmy(iy)+ddifpy(ny+1-iy) )/2.0
	  ddify(iy)  = ( ddify(iy)+ddify(ny+1-iy) )/2.0
	  meanpy(iy) = ( meanpy(iy)+meanmy(ny+1-iy) )/2.0
	  meanmy(iy) = ( meanmy(iy)+meanpy(ny+1-iy) )/2.0
  460   continue
	  y(iy0) = 0.0
	  ddifpy(iy0) = ( ddifpy(iy0)+ddifmy(iy0) )/2.0
	  ddifmy(iy0) = ( ddifpy(iy0)+ddifmy(iy0) )/2.0
	  ddify(iy0) =  ( ddifpy(iy0)+ddifmy(iy0) )/2.0
	  meanpy(iy0) = ( meanpy(iy0)+meanmy(iy0) )/2.0
	  meanmy(iy0) = ( meanpy(iy0)+meanmy(iy0) )/2.0
	do 470 iy = iy0+1, ny
	  y(iy)      = -y(ny+1-iy)
	  dify(iy)   = dify(ny+1-iy)
	  ddifpy(iy) = ddifmy(ny+1-iy)
	  ddifmy(iy) = ddifpy(ny+1-iy)
	  ddify(iy)  = ddify(ny+1-iy)
	  meanpy(iy) = meanmy(ny+1-iy)
	  meanmy(iy) = meanpy(ny+1-iy)
  470   continue
      end if
      do 500 iy = 1, ny
	y(iy) = yanf + y(iy)
  500 continue
c
c  grid boundary
c...............................................
      y(1)        = 2.*ymin - y(3)
      y(ny)       = 2.*ymax - y(ny-2)
      y(2)        = ymin
      y(ny1)     = ymax
      dify(1)     = dify(3)
      dify(ny)    = dify(ny-2)
      ddifmy(2)   = ( ddifpy(2) + ddifmy(2) ) / 2.0
      ddifmy(ny1)= ( ddifpy(ny1) + ddifmy(ny1) ) / 2.0
      ddifpy(2)   = ddifmy(2)
      ddifpy(ny1)= ddifmy(ny1)
      ddifmy(1)   = ddifpy(3)
      ddifpy(1)   = ddifmy(3)
      ddifmy(ny)  = ddifpy(ny-2)
      ddifpy(ny)  = ddifmy(ny-2)
      meanmy(1)   = meanpy(3)
      meanpy(1)   = meanmy(3)
      meanmy(ny)  = meanpy(ny-2)
      meanpy(ny)  = meanmy(ny-2)
c
      if (zentr(2))   then
	do 560 iy = 1,ny
	  iys = 1+ny-iy
	  if ( meanpy(iy) .ne. meanmy(iys) )
     +      write(26,31) iy,meanpy(iy),meanmy(iys),time
	  if ( ddifpy(iy) .ne. ddifmy(iys) )
     +      write(26,32) iy,ddifpy(iy),ddifmy(iys),time
  560   continue
      end if
   31 format(' meanpy(',i3,') =',f12.7,
     +       '   und',f12.7,'   time =',f9.4)
   32 format(' ddifpy(',i3,') =',f12.7,
     +       '   und',f12.7,'   time =',f9.4)

c.......................................................
c  z-coordinate
c      analogue  x
c.......................................................
      kz =  (zmax-zmin)/(nz-3)
      zanf = zmin
      iz0 = 2
      beta = pi/(zmax-zmin)
      if (zentr(3)) then
	zanf = (zmin+zmax)/2.0
	iz0 = (nz+1)/2
	beta = 2*pi/(zmax-zmin)
      end if
      a1 = 0.0
      a2 = 0.0
      a3 = 0.0
      a4 = 0.0
      if (.not. aequi(3)) then
	a1 = 24./17.*(eps(3)-kz)/kz/beta
	a2 =  6./17.*(kz-eps(3))/kz/beta
	a3 =  8./51.*(eps(3)-kz)/kz/beta
	a4 =  3./68.*(kz-eps(3))/kz/beta
      end if
c.......................................
      do 600 iz = 1,nz
         w      =  kz*(iz-iz0)
	 z(iz)  =  w + a1*sin(beta*w) + a2*sin(2*beta*w)
     +               + a3*sin(3.*beta*w) + a4*sin(4.*beta*w)
	 difz(iz) = 0.5/kz/(1. + a1*beta*cos(beta*w)
     +                    + 2.*a2*beta*cos(2.*beta*w)
     +                     + 3.*a3*beta*cos(3.*beta*w)
     +                      + 4.*a4*beta*cos(4.*beta*w) )
         w1 = w + 0.5*kz
         w2 = w - 0.5*kz
	 ddifpz(iz) = 2.*difz(iz)/kz/(1. + a1*beta*cos(beta*w1)
     +                              + 2.*a2*beta*cos(2.*beta*w1)
     +                               + 3.*a3*beta*cos(3.*beta*w1)
     +                                + 4.*a4*beta*cos(4.*beta*w1) )
	 ddifmz(iz) = 2.*difz(iz)/kz/(1. + a1*beta*cos(beta*w2)
     +                              + 2.*a2*beta*cos(2.*beta*w2)
     +                               + 3.*a3*beta*cos(3.*beta*w2)
     +                                + 4.*a4*beta*cos(4.*beta*w2) )
         ddifz(iz) = .5*(ddifpz(iz)+ddifmz(iz))
  600 continue
c
c  averages
c............................................
      do 650 iz = 2,nz1
	 meanmz(iz) = (z(iz+1)-z(iz))/(z(iz+1)-z(iz-1))
         meanpz(iz) = (z(iz)-z(iz-1))/(z(iz+1)-z(iz-1))
  650 continue
      meanmz(1)   = meanpz(3)
      meanpz(1)   = meanmz(3)
      meanmz(nz)  = meanpz(nz-2)
      meanpz(nz)  = meanmz(nz-2)
c
c centered nonuniform grid
c.......................................................
      if (zentr(3))   then
	do 660 iz = 1, iz0
	  z(iz)      = ( z(iz)-z(nz+1-iz) )/2.0
	  difz(iz)   = ( difz(iz)+difz(nz+1-iz) )/2.0
	  ddifpz(iz) = ( ddifpz(iz)+ddifmz(nz+1-iz) )/2.0
	  ddifmz(iz) = ( ddifmz(iz)+ddifpz(nz+1-iz) )/2.0
	  ddifz(iz)  = ( ddifz(iz)+ddifz(nz+1-iz) )/2.0
	  meanpz(iz) = ( meanpz(iz)+meanmz(nz+1-iz) )/2.0
	  meanmz(iz) = ( meanmz(iz)+meanpz(nz+1-iz) )/2.0
  660   continue
	  z(iz0) = 0.0
	  ddifpz(iz0) = ( ddifpz(iz0)+ddifmz(iz0) )/2.0
	  ddifmz(iz0) = ( ddifpz(iz0)+ddifmz(iz0) )/2.0
	  ddifz(iz0) =  ( ddifpz(iz0)+ddifmz(iz0) )/2.0
	  meanpz(iz0) = ( meanpz(iz0)+meanmz(iz0) )/2.0
	  meanmz(iz0) = ( meanpz(iz0)+meanmz(iz0) )/2.0
	do 670 iz = iz0+1, nz
	  z(iz)      = -z(nz+1-iz)
	  difz(iz)   = difz(nz+1-iz)
	  ddifpz(iz) = ddifmz(nz+1-iz)
	  ddifmz(iz) = ddifpz(nz+1-iz)
	  ddifz(iz)  = ddifz(nz+1-iz)
	  meanpz(iz) = meanmz(nz+1-iz)
	  meanmz(iz) = meanpz(nz+1-iz)
  670   continue
      end if
      do 700 iz = 1, nz
	z(iz) = zanf + z(iz)
  700 continue
c
c  grid boundary
c..............................................
      z(1)        = 2.*zmin - z(3)
      z(nz)       = 2.*zmax -z(nz-2)
      z(2)        = zmin
      z(nz1)     = zmax
      difz(1)     = difz(3)
      difz(nz)    = difz(nz-2)
      ddifmz(2)   = ( ddifpz(2) + ddifmz(2) ) / 2.0
      ddifmz(nz1)= ( ddifpz(nz1) + ddifmz(nz1) ) / 2.0
      ddifpz(2)   = ddifmz(2)
      ddifpz(nz1)= ddifmz(nz1)
      ddifmz(1)   = ddifpz(3)
      ddifpz(1)   = ddifmz(3)
      ddifmz(nz)  = ddifpz(nz-2)
      ddifpz(nz)  = ddifmz(nz-2)
      meanmz(1)   = meanpz(3)
      meanpz(1)   = meanmz(3)
      meanmz(nz)  = meanpz(nz-2)
      meanpz(nz)  = meanmz(nz-2)

      if (zentr(3)) then
	do 760 iz = 1,nz
	  if ( meanpz(iz) .ne. meanmz(nz+1-iz) )
     +      write(26,41) iz,meanpz(iz),meanmz(nz+1-iz),time
	  if ( ddifpz(iz) .ne. ddifmz(nz+1-iz) )
     +      write(26,42) iz,ddifpz(iz),ddifmz(nz+1-iz),time
  760   continue
      end if
   41 format(' meanpz(',i3,') =',f12.7,
     +       '   und',f12.7,'   time =',f9.4)
   42 format(' ddifpz(',i3,') =',f12.7,
     +       '   und',f12.7,'   time =',f9.4)
c
      return
      end
