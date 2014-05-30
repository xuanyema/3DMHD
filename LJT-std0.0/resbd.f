	subroutine resbd
c***********************************************************************
	include 'misflin'
c
	integer  ix,iy,iz,nx2,nx4,ny2,ny4,nz2,nz4,
     +           nxhalb,nyhalb,nzhalb
c ............................................................
c   computation of boundary values for the resist
c       f(i+1) = c1*f(i-1) + c2*f(i-3) + perio*f(n-3)
c                + lsym*f(sym)
c.............................................................
c     describtion in subroutine rdkoef
c     (line symmetry parameter lsym)
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
c--------------------------------------
c   x = xmin
c       (linesymmetry along z axis)
c--------------------------------------
      if ( perio(1) )  then
       do 110 iz = 2, nz1
       do 110 iy = 2, ny1
	 res(1,iy,iz) = res(3,iy,iz)
  110  continue
      else if ( lsym(1,1) ) then
       do 120 iz = 2, nz1
       do 120 iy = 2, ny1
	 res(1,iy,iz) = res(3,ny+1-iy,iz)
  120  continue
      else
       do 140 iz = 2, nz1
       do 140 iy = 2, ny1
	 res(1,iy,iz) = cres1(1,1)*res(3,iy,iz) + cres2(1,1)*res(5,iy,iz)
  140  continue
      end if
c--------------------------------------
c   x = xmax
c       (linesymmetry along z axis)
c--------------------------------------
      if ( perio(1) )  then
       do 210 iz = 2, nz1
       do 210 iy = 2, ny1
	res(nx,iy,iz) = res(3,iy,iz)
  210  continue
      else if ( lsym(2,1) )  then
       do 220 iz = 2, nz1
       do 220 iy = 2, ny1
	res(nx,iy,iz) =  res(nx2,ny+1-iy,iz)
  220  continue
      else
       do 240 iz = 2, nz1
       do 240 iy = 2, ny1
	res(nx,iy,iz) = cres1(2,1)*res(nx2,iy,iz) 
     +                  + cres2(2,1)*res(nx4,iy,iz)
  240  continue
      end if
c--------------------------------------
c   y = ymin
c       (linesymmetry along x axis)
c--------------------------------------
      if ( perio(2) )  then
       do 310 iz = 2, nz1
       do 310 ix = 1, nx
	res(ix,1,iz) = res(ix,ny2,iz)
  310  continue
      else if ( lsym(1,2) )  then
       do 320 iz = 2, nz1
       do 320 ix = 1, nx
	res(ix,1,iz) =  res(ix,3,nz+1-iz)
  320 continue
      else
       do 340 iz = 2, nz1
       do 340 ix = 1, nx
	res(ix,1,iz)   = cres1(1,2)*res(ix,3,iz)
     +                + cres2(1,2)*res(ix,5,iz)
  340  continue
      end if
c--------------------------------------
c   y = ymax
c       (linesymmetry along x axis)
c--------------------------------------
      if ( perio(2) )  then
       do 410 iz = 2, nz1
       do 410 ix = 1, nx
	res(ix,ny,iz) = res(ix,3,iz)
  410  continue
      else if ( lsym(2,2) )  then
       do 420 iz = 2, nz1
       do 420 ix = 1, nx
	res(ix,ny,iz) =  res(ix,ny2,nz+1-iz)
  420  continue
      else
       do 440 iz = 2, nz1
       do 440 ix = 1, nx
	res(ix,ny,iz) = cres1(2,2)*res(ix,ny2,iz)
     +                + cres2(2,2)*res(ix,ny4,iz)
  440  continue
      end if
c--------------------------------------
c   z = zmin
c--------------------------------------
      if ( perio(3) )  then
       do 510 iy = 1, ny
       do 510 ix = 1, nx
	res(ix,iy,1) = res(ix,iy,nz2)
  510  continue
      else if ( lsym(1,3) )  then
       do 520 iy = 1, ny
       do 520 ix = 1, nx
	res(ix,iy,1) =  res(ix,ny+1-iy,3)
c	res(ix,iy,1) =  res(nx+1-ix,iy,3)
  520  continue
      else
       do 540 iy = 1, ny
       do 540 ix = 1, nx
	res(ix,iy,1) = cres1(1,3)*res(ix,iy,3) + cres2(1,3)*res(ix,iy,5)
  540  continue
      end if
c--------------------------------------
c   z = zmax
c--------------------------------------
      if ( perio(3) )  then
       do 610 iy = 1, ny
       do 610 ix = 1, nx
	res(ix,iy,nz) = res(ix,iy,3)
  610  continue
      else if ( lsym(2,3) )  then
       do 620 iy = 1, ny
       do 620 ix = 1, nx
	res(ix,iy,nz) =  res(ix,ny+1-iy,nz2)
c	res(ix,iy,nz) =  res(nx+1-ix,iy,nz2)
  620  continue
      else
       do 640 iy = 1, ny
       do 640 ix = 1, nx
	res(ix,iy,nz) = cres1(2,3)*res(ix,iy,nz2) 
     +                  + cres2(2,3)*res(ix,iy,nz4)
  640  continue
      end if
c
      return
      end
