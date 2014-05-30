	subroutine bound
c***********************************************************************
c
	include 'misflin'
c
	integer  ix,iz,iy,nx2,nx4,ny2,ny4,nz2,nz4,
     +           nxhalb,nyhalb,nzhalb
	real     vxo,vzo
c ............................................................
c   computation of boundary values
c       f(i+1) = c1*f(i-1) + c2*f(i-3) + perio*f(n-3)
c                + lsym*f(sym)
c.............................................................
c     (line symmetry or periodicity determ by lsym or perio) 
c
      nxhalb = (nx+1)/2
      nyhalb = (ny+1)/2
      nzhalb = (nz+1)/2
      nx2 = nx-2
      nx4 = nx-4
      ny2 = ny-2
      ny4 = ny-4
      nz2 = nz-2
      nz4 = nz-4
c
c--------------------------------------
c   x = xmin
c       (linesymmetry along z axis)
c--------------------------------------
      if (1==0) then
      if ( perio(1) )  then

       do 110 iz = 2, nz1
       do 110 iy = 2, ny1
	 rho(1,iy,iz) = rho(nx2,iy,iz)
	 u(1,iy,iz)   = u(nx2,iy,iz)
	 res(1,iy,iz) = res(nx2,iy,iz)
	 bx(1,iy,iz)  = bx(nx2,iy,iz)
	 by(1,iy,iz)  = by(nx2,iy,iz)
	 bz(1,iy,iz)  = bz(nx2,iy,iz)
	 sx(1,iy,iz)  = sx(nx2,iy,iz)
	 sy(1,iy,iz)  = sy(nx2,iy,iz)
	 sz(1,iy,iz)  = sz(nx2,iy,iz)
  110  continue

      else if ( lsym(1,1) ) then

       do 120 iz = 2, nz1
       do 120 iy = 2, ny1
	 rho(1,iy,iz) = rho(3,ny+1-iy,iz)
	 u(1,iy,iz)   = u(3,ny+1-iy,iz)
	 res(1,iy,iz) = res(3,ny+1-iy,iz)
	 bx(1,iy,iz)  = bx(3,ny+1-iy,iz)
	 by(1,iy,iz)  = by(3,ny+1-iy,iz)
	 bz(1,iy,iz)  = -bz(3,ny+1-iy,iz)
	 sx(1,iy,iz)  = -sx(3,ny+1-iy,iz)
	 sy(1,iy,iz)  = -sy(3,ny+1-iy,iz)
	 sz(1,iy,iz)  = sz(3,ny+1-iy,iz)
  120  continue
       do 125 iz = 2, nz1
	 bz(2,nyhalb,iz)  = 0.0
	 sx(2,nyhalb,iz)  = 0.0
	 sy(2,nyhalb,iz)  = 0.0
  125  continue

      else

       do 140 iz = 2, nz1
       do 140 iy = 2, ny1
	 rho(1,iy,iz) = crho1(1,1)*rho(3,iy,iz) + crho2(1,1)*rho(5,iy,iz)
	 u(1,iy,iz)   = cu1(1,1)*u(3,iy,iz) + cu2(1,1)*u(5,iy,iz)
	 res(1,iy,iz) = cres1(1,1)*res(3,iy,iz) + cres2(1,1)*res(5,iy,iz)
	 bx(1,iy,iz)  = cbx1(1,1)*bx(3,iy,iz) + cbx2(1,1)*bx(5,iy,iz)
	 by(1,iy,iz)  = cby1(1,1)*by(3,iy,iz) + cby2(1,1)*by(5,iy,iz)
	 bz(1,iy,iz)  = cbz1(1,1)*bz(3,iy,iz) + cbz2(1,1)*bz(5,iy,iz)
	 sx(1,iy,iz)  = csx1(1,1)*sx(3,iy,iz) + csx2(1,1)*sx(5,iy,iz)
	 sy(1,iy,iz)  = csy1(1,1)*sy(3,iy,iz) + csy2(1,1)*sy(5,iy,iz)
	 sz(1,iy,iz)  = csz1(1,1)*sz(3,iy,iz) + csz2(1,1)*sz(5,iy,iz)
  140  continue

       if (ksx(1,1) .eq. -1.) then
	do 141 iz = 2, nz1
	do 141 iy = 2, ny1
  141     sx(2,iy,iz) = 0.
       end if
       if (ksy(1,1) .eq. -1.) then
	do 142 iz = 2, nz1
	do 142 iy = 2, ny1
  142     sy(2,iy,iz) = 0.
       end if
       if (ksz(1,1) .eq. -1.) then
	do 143 iz = 2, nz1
	do 143 iy = 2, ny1
  143     sz(2,iy,iz) = 0.
       end if
       if (kbx(1,1) .eq. -1.) then
	do 144 iz = 2, nz1
	do 144 iy = 2, ny1
  144	   bx(2,iy,iz) = 0.
       end if
       if (kby(1,1) .eq. -1.) then
	do 147 iz = 2, nz1
	do 147 iy = 2, ny1
  147     by(2,iy,iz) = 0.
       end if
       if (kbz(1,1) .eq. -1.) then
	do 148 iz = 2, nz1
	do 148 iy = 2, ny1
  148     bz(2,iy,iz) = 0.
       end if
      end if
      end if

      call boundxt


c--------------------------------------
c   x = xmax
c       (linesymmetry along z axis)
c--------------------------------------
      if ( perio(1) )  then

       do 210 iz = 2, nz1
       do 210 iy = 2, ny1
	rho(nx,iy,iz) = rho(3,iy,iz)
	u(nx,iy,iz)   = u(3,iy,iz)
	res(nx,iy,iz) = res(3,iy,iz)
	bx(nx,iy,iz)  = bx(3,iy,iz)
	by(nx,iy,iz)  = by(3,iy,iz)
	bz(nx,iy,iz)  = bz(3,iy,iz)
	sx(nx,iy,iz)  = sx(3,iy,iz)
	sy(nx,iy,iz)  = sy(3,iy,iz)
	sz(nx,iy,iz)  = sz(3,iy,iz)
  210  continue

      else if ( lsym(2,1) )  then

       do 220 iz = 2, nz1
       do 220 iy = 2, ny1
	rho(nx,iy,iz) =  rho(nx2,ny+1-iy,iz)
	u(nx,iy,iz)   =  u(nx2,ny+1-iy,iz)
	res(nx,iy,iz) =  res(nx2,ny+1-iy,iz)
	bx(nx,iy,iz)  =  bx(nx2,ny+1-iy,iz)
	by(nx,iy,iz)  =  by(nx2,ny+1-iy,iz)
	bz(nx,iy,iz)  = -bz(nx2,ny+1-iy,iz)
	sx(nx,iy,iz)  = -sx(nx2,ny+1-iy,iz)
	sy(nx,iy,iz)  = -sy(nx2,ny+1-iy,iz)
	sz(nx,iy,iz)  =  sz(nx2,ny+1-iy,iz)
  220  continue
       do 225 iz = 2, nz1
	 bz(nx1,nyhalb,iz)  = 0.0
	 sx(nx1,nyhalb,iz)  = 0.0
	 sy(nx1,nyhalb,iz)  = 0.0
  225  continue

      else

       do 240 iz = 2, nz1
       do 240 iy = 2, ny1
	rho(nx,iy,iz)= crho1(2,1)*rho(nx2,iy,iz) 
     +                 + crho2(2,1)*rho(nx4,iy,iz)
	u(nx,iy,iz)  = cu1(2,1)*u(nx2,iy,iz) + cu2(2,1)*u(nx4,iy,iz)
	res(nx,iy,iz)= cres1(2,1)*res(nx2,iy,iz) 
     +                 + cres2(2,1)*res(nx4,iy,iz)
	bx(nx,iy,iz) = bnxmax(1,iy,iz)
        bx(nx,iy,iz) = cbx1(2,1)*bx(nx2,iy,iz) + cbx2(2,1)*bx(nx4,iy,iz)
	by(nx,iy,iz) = cby1(2,1)*by(nx2,iy,iz) + cby2(2,1)*by(nx4,iy,iz)
	bz(nx,iy,iz) = cbz1(2,1)*bz(nx2,iy,iz) + cbz2(2,1)*bz(nx4,iy,iz)
	sx(nx,iy,iz) = csx1(2,1)*sx(nx2,iy,iz) + csx2(2,1)*sx(nx4,iy,iz)
	sy(nx,iy,iz) = csy1(2,1)*sy(nx2,iy,iz) + csy2(2,1)*sy(nx4,iy,iz)
	sz(nx,iy,iz) = csz1(2,1)*sz(nx2,iy,iz) + csz2(2,1)*sz(nx4,iy,iz)
  240  continue

       if (ksx(2,1) .eq. -1.) then
	do 241 iz = 2, nz1
	do 241 iy = 2, ny1
  241     sx(nx1,iy,iz) = 0.
       end if
       if (ksy(2,1) .eq. -1.) then
	do 242 iz = 2, nz1
	do 242 iy = 2, ny1
  242     sy(nx1,iy,iz) = 0.
       end if
       if (ksz(2,1) .eq. -1.) then
	do 243 iz = 2, nz1
	do 243 iy = 2, ny1
  243     sz(nx1,iy,iz) = 0.
       end if
       if (kbx(2,1) .eq. -1.) then
	do 250 iz = 2, nz1
	do 250 iy = 2, ny1
  250	   bx(nx1,iy,iz) = 0.
       end if
       if (kby(2,1) .eq. -1.) then
	do 255 iz = 2, nz1
	do 255 iy = 2, ny1
  255	   by(nx1,iy,iz) = 0.
       end if
       if (kbz(2,1) .eq. -1.) then
	do 260 iz = 2, nz1
	do 260 iy = 2, ny1
  260	   bz(nx1,iy,iz) = 0.
       end if


      end if

c--------------------------------------
c   y = ymin
c       (linesymmetry along x axis)
c--------------------------------------
      if ( perio(2) )  then

       do 310 iz = 2, nz1
       do 310 ix = 1, nx
	rho(ix,1,iz) = rho(ix,ny2,iz)
	u(ix,1,iz)   = u(ix,ny2,iz)
	res(ix,1,iz) = res(ix,ny2,iz)
	bx(ix,1,iz)  = bx(ix,ny2,iz)
	by(ix,1,iz)  = by(ix,ny2,iz)
	bz(ix,1,iz)  = bz(ix,ny2,iz)
	sx(ix,1,iz)  = sx(ix,ny2,iz)
	sy(ix,1,iz)  = sy(ix,ny2,iz)
	sz(ix,1,iz)  = sz(ix,ny2,iz)
  310  continue

      else if ( lsym(1,2) )  then

       do 320 iz = 2, nz1
       do 320 ix = 1, nx
	rho(ix,1,iz) =  rho(ix,3,nz+1-iz)
	u(ix,1,iz)   =  u(ix,3,nz+1-iz)
	res(ix,1,iz) =  res(ix,3,nz+1-iz)
	bx(ix,1,iz)  = -bx(ix,3,nz+1-iz)
	by(ix,1,iz)  =  by(ix,3,nz+1-iz)
	bz(ix,1,iz)  =  bz(ix,3,nz+1-iz)
	sx(ix,1,iz)  =  sx(ix,3,nz+1-iz)
	sy(ix,1,iz)  = -sy(ix,3,nz+1-iz)
	sz(ix,1,iz)  = -sz(ix,3,nz+1-iz)
  320 continue
       do 325 ix = 1, nx
	bx(ix,2,nzhalb)  = 0.0
	sy(ix,2,nzhalb)  = 0.0
	sz(ix,2,nzhalb)  = 0.0
  325  continue

      else

       do 340 iz = 2, nz1
       do 340 ix = 1, nx
	rho(ix,1,iz) = crho1(1,2)*rho(ix,3,iz)
     +                + crho2(1,2)*rho(ix,5,iz)
	u(ix,1,iz)   = cu1(1,2)*u(ix,3,iz) + cu2(1,2)*u(ix,5,iz)
	res(ix,1,iz)   = cres1(1,2)*res(ix,3,iz)
     +                + cres2(1,2)*res(ix,5,iz)
	bx(ix,1,iz)  = cbx1(1,2)*bx(ix,3,iz)+cbx2(1,2)*bx(ix,5,iz)
	by(ix,1,iz)  = cby1(1,2)*by(ix,3,iz)+cby2(1,2)*by(ix,5,iz)
	bz(ix,1,iz)  = cbz1(1,2)*bz(ix,3,iz)+cbz2(1,2)*bz(ix,5,iz)
	sx(ix,1,iz)  = csx1(1,2)*sx(ix,3,iz)+csx2(1,2)*sx(ix,5,iz)
 	sy(ix,1,iz)  = csy1(1,2)*sy(ix,3,iz)+csy2(1,2)*sy(ix,5,iz)
	sz(ix,1,iz)  = csz1(1,2)*sz(ix,3,iz)+csz2(1,2)*sz(ix,5,iz)
  340  continue

       if (ksx(1,2) .eq. -1.) then
	do 341 iz = 2, nz1
	do 341 ix = 1, nx
  341     sx(ix,2,iz) = 0.
       end if
       if (ksy(1,2) .eq. -1.) then
	do 342 iz = 2, nz1
	do 342 ix = 1, nx
  342     sy(ix,2,iz) = 0.
       end if
       if (ksz(1,2) .eq. -1.) then
	do 343 iz = 2, nz1
	do 343 ix = 1, nx
  343     sz(ix,2,iz) = 0.
       end if
       if (kbx(1,2) .eq. -1.) then
	do 346 iz = 2, nz1
	do 346 ix = 1, nx
  346     bx(ix,2,iz) = 0.
       end if
       if (kby(1,2) .eq. -1.) then
	do 347 iz = 2, nz1
	do 347 ix = 1, nx
  347	   by(ix,2,iz) = 0.
       end if
       if (kbz(1,2) .eq. -1.) then
	do 349 iz = 2, nz1
	do 349 ix = 1, nx
  349	   bz(ix,2,iz) = 0.
       end if


      end if

c--------------------------------------
c   y = ymax
c       (linesymmetry along x axis)
c--------------------------------------
      if ( perio(2) )  then

       do 410 iz = 2, nz1
       do 410 ix = 1, nx
	rho(ix,ny,iz) = rho(ix,3,iz)
	u(ix,ny,iz)   = u(ix,3,iz)
	res(ix,ny,iz) = res(ix,3,iz)
	bx(ix,ny,iz)  = bx(ix,3,iz)
	by(ix,ny,iz)  = by(ix,3,iz)
	bz(ix,ny,iz)  = bz(ix,3,iz)
	sx(ix,ny,iz)  = sx(ix,3,iz)
	sy(ix,ny,iz)  = sy(ix,3,iz)
	sz(ix,ny,iz)  = sz(ix,3,iz)
  410  continue

      else if ( lsym(2,2) )  then

       do 420 iz = 2, nz1
       do 420 ix = 1, nx
	rho(ix,ny,iz) =  rho(ix,ny2,nz+1-iz)
	u(ix,ny,iz)   =  u(ix,ny2,nz+1-iz)
	res(ix,ny,iz) =  res(ix,ny2,nz+1-iz)
	bx(ix,ny,iz)  = -bx(ix,ny2,nz+1-iz)
	by(ix,ny,iz)  =  by(ix,ny2,nz+1-iz)
	bz(ix,ny,iz)  =  bz(ix,ny2,nz+1-iz)
	sx(ix,ny,iz)  =  sx(ix,ny2,nz+1-iz)
 	sy(ix,ny,iz)  = -sy(ix,ny2,nz+1-iz)
	sz(ix,ny,iz)  = -sz(ix,ny2,nz+1-iz)
  420  continue
       do 425 ix = 1, nx
	bx(ix,ny1,nzhalb)  = 0.0
	sy(ix,ny1,nzhalb)  = 0.0
	sz(ix,ny1,nzhalb)  = 0.0
  425  continue

      else

       do 440 iz = 2, nz1
       do 440 ix = 1, nx
	rho(ix,ny,iz) = crho1(2,2)*rho(ix,ny2,iz)
     +                + crho2(2,2)*rho(ix,ny4,iz)
	u(ix,ny,iz)   = cu1(2,2)*u(ix,ny2,iz)+cu2(2,2)*u(ix,ny4,iz)
	res(ix,ny,iz) = cres1(2,2)*res(ix,ny2,iz)
     +                + cres2(2,2)*res(ix,ny4,iz)
	bx(ix,ny,iz)  = cbx1(2,2)*bx(ix,ny2,iz)+cbx2(2,2)*bx(ix,ny4,iz)
	by(ix,ny,iz)  = cby1(2,2)*by(ix,ny2,iz)+cby2(2,2)*by(ix,ny4,iz)
	bz(ix,ny,iz)  = cbz1(2,2)*bz(ix,ny2,iz)+cbz2(2,2)*bz(ix,ny4,iz)
	sx(ix,ny,iz)  = csx1(2,2)*sx(ix,ny2,iz)+csx2(2,2)*sx(ix,ny4,iz)
	sy(ix,ny,iz)  = csy1(2,2)*sy(ix,ny2,iz)+csy2(2,2)*sy(ix,ny4,iz)
	sz(ix,ny,iz)  = csz1(2,2)*sz(ix,ny2,iz)+csz2(2,2)*sz(ix,ny4,iz)
  440  continue

       if (ksx(2,2) .eq. -1.) then
	do 441 iz = 2, nz1
	do 441 ix = 1, nx
  441     sx(ix,ny1,iz) = 0.
       end if
       if (ksy(2,2) .eq. -1.) then
	do 442 iz = 2, nz1
	do 442 ix = 1, nx
  442     sy(ix,ny1,iz) = 0.
       end if
       if (ksz(2,2) .eq. -1.) then
	do 443 iz = 2, nz1
	do 443 ix = 1, nx
  443     sz(ix,ny1,iz) = 0.
       end if
       if (kbx(2,2) .eq. -1.) then
	do 446 iz = 2, nz1
	do 446 ix = 1, nx
  446     bx(ix,ny1,iz) = 0.
       end if
       if (kby(2,2) .eq. -1.) then
	do 447 iz = 2, nz1
	do 447 ix = 1, nx
  447     by(ix,ny1,iz) = 0.
c      else
c       do 448 iz = 2, nz1
c       do 448 ix = 1, nx
c 448	  by(ix,ny1,iz) = bnymax(ix,2,iz)
       end if
       if (kbz(2,2) .eq. -1.) then
	do 449 iz = 2, nz1
	do 449 ix = 1, nx
  449	   bz(ix,ny1,iz) = 0.
       end if


      end if

c--------------------------------------
c   z = zmin
c--------------------------------------
      if ( perio(3) )  then

       do 510 iy = 1, ny
       do 510 ix = 1, nx
	rho(ix,iy,1) = rho(ix,iy,nz2)
	u(ix,iy,1)   = u(ix,iy,nz2)
	res(ix,iy,1) = res(ix,iy,nz2)
	bx(ix,iy,1)  = bx(ix,iy,nz2)
	by(ix,iy,1)  = by(ix,iy,nz2)
	bz(ix,iy,1)  = bz(ix,iy,nz2)
	sx(ix,iy,1)  = sx(ix,iy,nz2)
	sy(ix,iy,1)  = sy(ix,iy,nz2)
	sz(ix,iy,1)  = sz(ix,iy,nz2)
  510  continue

      else if ( lsym(1,3) )  then

c for magnetopause computation linesymmetry along x
       do 520 iy = 1, ny
       do 520 ix = 1, nx
	rho(ix,iy,1) =  rho(ix,ny+1-iy,3)
	u(ix,iy,1)   =  u(ix,ny+1-iy,3)
	res(ix,iy,1) =  res(ix,ny+1-iy,3)
	bx(ix,iy,1)  = -bx(ix,ny+1-iy,3)
	by(ix,iy,1)  =  by(ix,ny+1-iy,3)
	bz(ix,iy,1)  =  bz(ix,ny+1-iy,3)
	sx(ix,iy,1)  =  sx(ix,ny+1-iy,3)
	sy(ix,iy,1)  = -sy(ix,ny+1-iy,3)
	sz(ix,iy,1)  = -sz(ix,ny+1-iy,3)
  520  continue
       do 525 ix = 1,nx
	bx(ix,nyhalb,2)  = 0.0
	sy(ix,nyhalb,2)  = 0.0
	sz(ix,nyhalb,2)  = 0.0
  525  continue

c for arcs ( y(pro)->z, z(pro)->-y ) linesym along y(pro):
c       do 520 iy = 1, ny
c       do 520 ix = 1, nx
c	rho(ix,iy,1) =  rho(nx+1-ix,iy,3)
c	u(ix,iy,1)   =  u(nx+1-ix,iy,3)
c	res(ix,iy,1) =  res(nx+1-ix,iy,3)
c	bx(ix,iy,1)  = -bx(nx+1-ix,iy,3)
c	by(ix,iy,1)  =  by(nx+1-ix,iy,3)
c	bz(ix,iy,1)  = -bz(nx+1-ix,iy,3)
c	sx(ix,iy,1)  = -sx(nx+1-ix,iy,3)
c	sy(ix,iy,1)  =  sy(nx+1-ix,iy,3)
c	sz(ix,iy,1)  = -sz(nx+1-ix,iy,3)
c  520  continue
c       do 525 iy = 1,ny
c	bx(nxhalb,iy,2)  = 0.0
c	bz(nxhalb,iy,2)  = 0.0
c	sx(nxhalb,iy,2)  = 0.0
c	sz(nxhalb,iy,2)  = 0.0
c  525  continue

      else

       do 540 iy = 1, ny
       do 540 ix = 1, nx
	rho(ix,iy,1)= crho1(1,3)*rho(ix,iy,3) + crho2(1,3)*rho(ix,iy,5)
	u(ix,iy,1)  = cu1(1,3)*u(ix,iy,3) + cu2(1,3)*u(ix,iy,5)
	res(ix,iy,1)= cres1(1,3)*res(ix,iy,3) + cres2(1,3)*res(ix,iy,5)
	bx(ix,iy,1) = cbx1(1,3)*bx(ix,iy,3) + cbx2(1,3)*bx(ix,iy,5)
	by(ix,iy,1) = cby1(1,3)*by(ix,iy,3) + cby2(1,3)*by(ix,iy,5)
	bz(ix,iy,1) = cbz1(1,3)*bz(ix,iy,3) + cbz2(1,3)*bz(ix,iy,5)
	sx(ix,iy,1) = csx1(1,3)*sx(ix,iy,3) + csx2(1,3)*sx(ix,iy,5)
	sy(ix,iy,1) = csy1(1,3)*sy(ix,iy,3) + csy2(1,3)*sy(ix,iy,5)
	sz(ix,iy,1) = csz1(1,3)*sz(ix,iy,3) + csz2(1,3)*sz(ix,iy,5)
  540  continue

       if (ksx(1,3) .eq. -1.) then
	do 541 iy = 1, ny
	do 541 ix = 1, nx
  541     sx(ix,iy,2) = 0.
       end if
       if (ksy(1,3) .eq. -1.) then
	do 542 iy = 1, ny
	do 542 ix = 1, nx
  542     sy(ix,iy,2) = 0.
       end if
       if (ksz(1,3) .eq. -1.) then
	do 543 iy = 1, ny
	do 543 ix = 1, nx
  543     sz(ix,iy,2) = 0.
       end if
       if (kbx(1,3) .eq. -1.) then
	do 546 iy = 1, ny
	do 546 ix = 1, nx
  546     bx(ix,iy,2) = 0.
       end if
       if (kby(1,3) .eq. -1.) then
	do 547 iy = 1, ny
	do 547 ix = 1, nx
  547     by(ix,iy,2) = 0.
       end if
       if (kbz(1,3) .eq. -1.) then
	do 548 iy = 1, ny
	do 548 ix = 1, nx
  548     bz(ix,iy,2) = 0.
       end if


      end if

c--------------------------------------
c   z = zmax
c--------------------------------------
      if ( perio(3) )  then

       do 610 iy = 1, ny
       do 610 ix = 1, nx
	rho(ix,iy,nz) = rho(ix,iy,3)
	u(ix,iy,nz)   = u(ix,iy,3)
	res(ix,iy,nz) = res(ix,iy,3)
	bx(ix,iy,nz)  = bx(ix,iy,3)
	by(ix,iy,nz)  = by(ix,iy,3)
	bz(ix,iy,nz)  = bz(ix,iy,3)
	sx(ix,iy,nz)  = sx(ix,iy,3)
	sy(ix,iy,nz)  = sy(ix,iy,3)
	sz(ix,iy,nz)  = sz(ix,iy,3)
  610  continue

      else if ( lsym(2,3) )  then

c for magnetopause computation linesymmetry along x
       do 620 iy = 1, ny
       do 620 ix = 1, nx
	rho(ix,iy,nz) =  rho(ix,ny+1-iy,nz2)
	u(ix,iy,nz)   =  u(ix,ny+1-iy,nz2)
	res(ix,iy,nz) =  res(ix,ny+1-iy,nz2)
	bx(ix,iy,nz)  = -bx(ix,ny+1-iy,nz2)
	by(ix,iy,nz)  =  by(ix,ny+1-iy,nz2)
	bz(ix,iy,nz)  =  bz(ix,ny+1-iy,nz2)
	sx(ix,iy,nz)  =  sx(ix,ny+1-iy,nz2)
	sy(ix,iy,nz)  = -sy(ix,ny+1-iy,nz2)
	sz(ix,iy,nz)  = -sz(ix,ny+1-iy,nz2)
  620  continue
       do 625 ix = 1,nx
	bx(ix,nyhalb,nz1)  = 0.0
	sy(ix,nyhalb,nz1)  = 0.0
	sz(ix,nyhalb,nz1)  = 0.0
  625  continue

c       do 620 iy = 1, ny
c       do 620 ix = 1, nx
c	rho(ix,iy,nz) =  rho(nx+1-ix,iy,nz2)
c	u(ix,iy,nz)   =  u(nx+1-ix,iy,nz2)
c	res(ix,iy,nz) =  res(nx+1-ix,iy,nz2)
c	bx(ix,iy,nz)  = -bx(nx+1-ix,iy,nz2)
c	by(ix,iy,nz)  =  by(nx+1-ix,iy,nz2)
c	bz(ix,iy,nz)  = -bz(nx+1-ix,iy,nz2)
c	sx(ix,iy,nz)  = -sx(nx+1-ix,iy,nz2)
c	sy(ix,iy,nz)  =  sy(nx+1-ix,iy,nz2)
c	sz(ix,iy,nz)  = -sz(nx+1-ix,iy,nz2)
c  620  continue
c       do 625 iy = 1,ny
c	bx(nxhalb,iy,nz1)  = 0.0
c	bz(nxhalb,iy,nz1)  = 0.0
c	sx(nxhalb,iy,nz1)  = 0.0
c	sz(nxhalb,iy,nz1)  = 0.0
c  625  continue

      else

       do 640 iy = 1, ny
       do 640 ix = 1, nx
	rho(ix,iy,nz)= crho1(2,3)*rho(ix,iy,nz2) 
     +                 + crho2(2,3)*rho(ix,iy,nz4)
	u(ix,iy,nz)  = cu1(2,3)*u(ix,iy,nz2) + cu2(2,3)*u(ix,iy,nz4)
	res(ix,iy,nz)= cres1(2,3)*res(ix,iy,nz2) 
     +                 + cres2(2,3)*res(ix,iy,nz4)
	bx(ix,iy,nz) = cbx1(2,3)*bx(ix,iy,nz2) + cbx2(2,3)*bx(ix,iy,nz4)
	by(ix,iy,nz) = cby1(2,3)*by(ix,iy,nz2) + cby2(2,3)*by(ix,iy,nz4)
	bz(ix,iy,nz) = cbz1(2,3)*bz(ix,iy,nz2) + cbz2(2,3)*bz(ix,iy,nz4)
	sx(ix,iy,nz) = csx1(2,3)*sx(ix,iy,nz2) + csx2(2,3)*sx(ix,iy,nz4)
	sy(ix,iy,nz) = csy1(2,3)*sy(ix,iy,nz2) + csy2(2,3)*sy(ix,iy,nz4)
	sz(ix,iy,nz) = csz1(2,3)*sz(ix,iy,nz2) + csz2(2,3)*sz(ix,iy,nz4)
  640  continue

       if (ksx(2,3) .eq. -1.) then
	do 641 iy = 1, ny
	do 641 ix = 1, nx
  641     sx(ix,iy,nz1) = 0.
       end if
       if (ksy(2,3) .eq. -1.) then
	do 642 iy = 1, ny
	do 642 ix = 1, nx
  642     sy(ix,iy,nz1) = 0.
       end if
       if (ksz(2,3) .eq. -1.) then
	do 643 iy = 1, ny
	do 643 ix = 1, nx
  643     sz(ix,iy,nz1) = 0.
       end if
       if (kbx(2,3) .eq. -1.) then
	do 646 iy = 1, ny
	do 646 ix = 1, nx
  646     bx(ix,iy,nz1) = 0.
       end if
       if (kby(2,3) .eq. -1.) then
	do 647 iy = 1, ny
	do 647 ix = 1, nx
  647     by(ix,iy,nz1) = 0.
       end if
       if (kbz(2,3) .eq. -1.) then
	do 648 iy = 1, ny
	do 648 ix = 1, nx
  648     bz(ix,iy,nz1) = 0.
       end if

      end if
c------------------------------------------------
c  div b=0 ! attention at edges
c------------------------------------------------
c
c     div b for bx 
c
c      a) bound cond at xmin,xmax
c      b) edges x,y = const
c      c) corners at zmin,zmax
c------------------------------------------------
c
c   x = xmin bzw xmax
c................................................
      if ( .not.(perio(1) .or. lsym(1,1)) )  then
       do 810 iz = 2, nz1
       do 810 iy = 2, ny1
	  bx(1,iy,iz)  = bx(3,iy,iz) + 1./difx(2) * (
     +                + dify(iy)*(by(3,iy+1,iz)-by(3,iy-1,iz))
     +                + difz(iz)*(bz(3,iy,iz+1)-bz(3,iy,iz-1)) )
  810  continue
      end if

c
      if ( .not.(perio(1) .or. lsym(2,1)) )  then
       do 820 iz = 2, nz1
       do 820 iy = 2, ny1
 	  bx(nx,iy,iz) = bx(nx2,iy,iz) - 1./difx(nx1) * (
     +                + dify(iy)*(by(nx1,iy+1,iz)-by(nx1,iy-1,iz))
     +                + difz(iz)*(bz(nx1,iy,iz+1)-bz(nx1,iy,iz-1)) )
  820  continue
      end if
c
c   y = ymin bzw ymax
c................................................
c      if ( .not.(perio(2) .or. lsym(1,2)) )  then
c       do 830 iz = 2, nz1
c       do 830 ix = 2, nx1
c	  by(ix,1,iz) = by(ix,3,iz) + 1./dify(2) * (
c     +                + difz(iz)*(bz(ix,2,iz+1)-bz(ix,2,iz-1))
c     +                + difx(ix)*(bx(ix+1,2,iz)-bx(ix-1,2,iz)) )
c  830  continue
c      end if
c
c      if ( .not.(perio(2) .or. lsym(2,2)) )  then
c       do 840 iz = 2, nz1
c       do 840 ix = 2, nx1
c	  by(ix,ny,iz) = by(ix,ny2,iz) - 1./dify(ny1) * (
c     +                + difz(iz)*(bz(ix,ny1,iz+1)-bz(ix,ny1,iz-1))
c     +                + difx(ix)*(bx(ix+1,ny1,iz)-bx(ix-1,ny1,iz)) )
c  840  continue
c      end if
c
c    z = zmin bzw zmax
c................................................
      if ( .not.(perio(3) .or. lsym(1,3)) )  then
       do 860 iy = 2, ny1
       do 860 ix = 2, nx1
	 bz(ix,iy,1)  = bz(ix,iy,3) + 1./difz(2) * (
     +               + difx(ix)*(bx(ix+1,iy,2)-bx(ix-1,iy,2))
     +               + dify(iy)*(by(ix,iy+1,2)-by(ix,iy-1,2)) )
  860  continue
      end if

      if ( .not.(perio(3) .or. lsym(2,3)) )  then
       do 862 iy = 2, ny1
       do 862 ix = 2, nx1
	 bz(ix,iy,nz)  = bz(ix,iy,nz2) - 1./difz(nz1) * (
     +                + difx(ix)*(bx(ix+1,iy,nz1)-bx(ix-1,iy,nz1))
     +                + dify(iy)*(by(ix,iy+1,nz1)-by(ix,iy-1,nz1)) )
  862  continue
      end if
c
      return
      end
