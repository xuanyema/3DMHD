	subroutine bdcoef
c***********************************************************************
	include 'misflin'
c
	integer  i,j,nx3,ny3,nz3
c ............................................................
c     bdcoef determiines coefficients for the boundary 
c     conditions in the following manner 
c ............................................................
c     f(n) = k*f(n-2) + alfa*f'(n-3)*del
c     mit
c        del = 1/dif(n-1)
c        f'(n-2)  = dif(n-3)*[f(n-2) - f(n-4)]
c   => f(n) = ( k + alfa*dif(n-3)/dif(n-1) )*f(n-2)
c                     - alfa*dif(n-3)/dif(n-1)*f(n-4)
c   und
c      f(1) = ( k + alfa*dif(4)/dif(2) )*f(3)
c                     - alfa*dif(4)/dif(2)*f(5)
c.............................................
c
c     antisymmetry      k = -1  ,  alfa = 0
c     symmetry          k = 1   ,  alfa = 0
c     extrapolation      k = 1   ,  alfa = 1
c
c     b from div b = 0 
c.............................................
c     boundcoeffizienten:
c       1. index:  1=min, 2=max
c       2. index:  1=x, 2=y, 3=z
c
      nx3 = nx-3
      ny3 = ny-3
      nz3 = nz-3
c
c..................
c     x = xmin
c..................
      crho1(1,1) = 1. + arho(1,1)*difx(4)/difx(2)
      crho2(1,1) =  - arho(1,1)*difx(4)/difx(2)
      cu1(1,1)   = 1. + au(1,1)*difx(4)/difx(2)
      cu2(1,1)   =  - au(1,1)*difx(4)/difx(2)
      cres1(1,1)   = 1. + ares(1,1)*difx(4)/difx(2)
      cres2(1,1)   =  - ares(1,1)*difx(4)/difx(2)
      cbx1(1,1)  = kbx(1,1) + abx(1,1)*difx(4)/difx(2)
      cbx2(1,1)  =  - abx(1,1)*difx(4)/difx(2)
      cby1(1,1)  = kby(1,1) + aby(1,1)*difx(4)/difx(2)
      cby2(1,1)  =  - aby(1,1)*difx(4)/difx(2)
      cbz1(1,1)  = kbz(1,1) + abz(1,1)*difx(4)/difx(2)
      cbz2(1,1)  =  - abz(1,1)*difx(4)/difx(2)
      csx1(1,1)  = ksx(1,1) + asx(1,1)*difx(4)/difx(2)
      csx2(1,1)  =  - asx(1,1)*difx(4)/difx(2)
      csy1(1,1)  = ksy(1,1) + asy(1,1)*difx(4)/difx(2)
      csy2(1,1)  =  - asy(1,1)*difx(4)/difx(2)
      csz1(1,1)  = ksz(1,1) + asz(1,1)*difx(4)/difx(2)
      csz2(1,1)  =  - asz(1,1)*difx(4)/difx(2)
c..................
c     x = xmax
c..................
      crho1(2,1) = 1. + arho(2,1)*difx(nx3)/difx(nx1)
      crho2(2,1) =  - arho(2,1)*difx(nx3)/difx(nx1)
      cu1(2,1)   = 1. + au(2,1)*difx(nx3)/difx(nx1)
      cu2(2,1)   =  - au(2,1)*difx(nx3)/difx(nx1)
      cres1(2,1)   = 1. + ares(2,1)*difx(nx3)/difx(nx1)
      cres2(2,1)   =  - ares(2,1)*difx(nx3)/difx(nx1)
      cbx1(2,1)  = kbx(2,1) + abx(2,1)*difx(nx3)/difx(nx1)
      cbx2(2,1)  =  - abx(2,1)*difx(nx3)/difx(nx1)
      cby1(2,1)  = kby(2,1) + aby(2,1)*difx(nx3)/difx(nx1)
      cby2(2,1)  =  - aby(2,1)*difx(nx3)/difx(nx1)
      cbz1(2,1)  = kbz(2,1) + abz(2,1)*difx(nx3)/difx(nx1)
      cbz2(2,1)  =  - abz(2,1)*difx(nx3)/difx(nx1)
      csx1(2,1)  = ksx(2,1) + asx(2,1)*difx(nx3)/difx(nx1)
      csx2(2,1)  =  - asx(2,1)*difx(nx3)/difx(nx1)
      csy1(2,1)  = ksy(2,1) + asy(2,1)*difx(nx3)/difx(nx1)
      csy2(2,1)  =  - asy(2,1)*difx(nx3)/difx(nx1)
      csz1(2,1)  = ksz(2,1) + asz(2,1)*difx(nx3)/difx(nx1)
      csz2(2,1)  =  - asz(2,1)*difx(nx3)/difx(nx1)
c..................
c     y = ymin
c..................
      crho1(1,2) = 1. + arho(1,2)*dify(4)/dify(2)
      crho2(1,2) =  - arho(1,2)*dify(4)/dify(2)
      cu1(1,2)   = 1. + au(1,2)*dify(4)/dify(2)
      cu2(1,2)   =  - au(1,2)*dify(4)/dify(2)
      cres1(1,2)   = 1. + ares(1,2)*dify(4)/dify(2)
      cres2(1,2)   =  - ares(1,2)*dify(4)/dify(2)
      cbx1(1,2)  = kbx(1,2) + abx(1,2)*dify(4)/dify(2)
      cbx2(1,2)  =  - abx(1,2)*dify(4)/dify(2)
      cby1(1,2)  = kby(1,2) + aby(1,2)*dify(4)/dify(2)
      cby2(1,2)  =  - aby(1,2)*dify(4)/dify(2)
      cbz1(1,2)  = kbz(1,2) + abz(1,2)*dify(4)/dify(2)
      cbz2(1,2)  =  - abz(1,2)*dify(4)/dify(2)
      csx1(1,2)  = ksx(1,2) + asx(1,2)*dify(4)/dify(2)
      csx2(1,2)  =  - asx(1,2)*dify(4)/dify(2)
      csy1(1,2)  = ksy(1,2) + asy(1,2)*dify(4)/dify(2)
      csy2(1,2)  =  - asy(1,2)*dify(4)/dify(2)
      csz1(1,2)  = ksz(1,2) + asz(1,2)*dify(4)/dify(2)
      csz2(1,2)  =  - asz(1,2)*dify(4)/dify(2)
c..................
c     y = ymax
c..................
      crho1(2,2) = 1. + arho(2,2)*dify(ny3)/dify(ny1)
      crho2(2,2) =  - arho(2,2)*dify(ny3)/dify(ny1)
      cu1(2,2)   = 1. + au(2,2)*dify(ny3)/dify(ny1)
      cu2(2,2)   =  - au(2,2)*dify(ny3)/dify(ny1)
      cres1(2,2)   = 1. + ares(2,2)*dify(ny3)/dify(ny1)
      cres2(2,2)   =  - ares(2,2)*dify(ny3)/dify(ny1)
      cbx1(2,2)  = kbx(2,2) + abx(2,2)*dify(ny3)/dify(ny1)
      cbx2(2,2)  =  - abx(2,2)*dify(ny3)/dify(ny1)
      cby1(2,2)  = kby(2,2) + aby(2,2)*dify(ny3)/dify(ny1)
      cby2(2,2)  =  - aby(2,2)*dify(ny3)/dify(ny1)
      cbz1(2,2)  = kbz(2,2) + abz(2,2)*dify(ny3)/dify(ny1)
      cbz2(2,2)  =  - abz(2,2)*dify(ny3)/dify(ny1)
      csx1(2,2)  = ksx(2,2) + asx(2,2)*dify(ny3)/dify(ny1)
      csx2(2,2)  =  - asx(2,2)*dify(ny3)/dify(ny1)
      csy1(2,2)  = ksy(2,2) + asy(2,2)*dify(ny3)/dify(ny1)
      csy2(2,2)  =  - asy(2,2)*dify(ny3)/dify(ny1)
      csz1(2,2)  = ksz(2,2) + asz(2,2)*dify(ny3)/dify(ny1)
      csz2(2,2)  =  - asz(2,2)*dify(ny3)/dify(ny1)
c..................
c     z = zmin
c..................
      crho1(1,3) = 1. + arho(1,3)*difz(4)/difz(2)
      crho2(1,3) =  - arho(1,3)*difz(4)/difz(2)
      cu1(1,3)   = 1. + au(1,3)*difz(4)/difz(2)
      cu2(1,3)   =  - au(1,3)*difz(4)/difz(2)
      cres1(1,3)   = 1. + ares(1,3)*difz(4)/difz(2)
      cres2(1,3)   =  - ares(1,3)*difz(4)/difz(2)
      cbx1(1,3)  = kbx(1,3) + abx(1,3)*difz(4)/difz(2)
      cbx2(1,3)  =  - abx(1,3)*difz(4)/difz(2)
      cby1(1,3)  = kby(1,3) + aby(1,3)*difz(4)/difz(2)
      cby2(1,3)  =  - aby(1,3)*difz(4)/difz(2)
      cbz1(1,3)  = kbz(1,3) + abz(1,3)*difz(4)/difz(2)
      cbz2(1,3)  =  - abz(1,3)*difz(4)/difz(2)
      csx1(1,3)  = ksx(1,3) + asx(1,3)*difz(4)/difz(2)
      csx2(1,3)  =  - asx(1,3)*difz(4)/difz(2)
      csy1(1,3)  = ksy(1,3) + asy(1,3)*difz(4)/difz(2)
      csy2(1,3)  =  - asy(1,3)*difz(4)/difz(2)
      csz1(1,3)  = ksz(1,3) + asz(1,3)*difz(4)/difz(2)
      csz2(1,3)  =  - asz(1,3)*difz(4)/difz(2)
c..................
c     z = zmax
c..................
      crho1(2,3) = 1. + arho(2,3)*difz(nz3)/difz(nz1)
      crho2(2,3) =  - arho(2,3)*difz(nz3)/difz(nz1)
      cu1(2,3)   = 1. + au(2,3)*difz(nz3)/difz(nz1)
      cu2(2,3)   =  - au(2,3)*difz(nz3)/difz(nz1)
      cres1(2,3)   = 1. + ares(2,3)*difz(nz3)/difz(nz1)
      cres2(2,3)   =  - ares(2,3)*difz(nz3)/difz(nz1)
      cbx1(2,3)  = kbx(2,3) + abx(2,3)*difz(nz3)/difz(nz1)
      cbx2(2,3)  =  - abx(2,3)*difz(nz3)/difz(nz1)
      cby1(2,3)  = kby(2,3) + aby(2,3)*difz(nz3)/difz(nz1)
      cby2(2,3)  =  - aby(2,3)*difz(nz3)/difz(nz1)
      cbz1(2,3)  = kbz(2,3) + abz(2,3)*difz(nz3)/difz(nz1)
      cbz2(2,3)  =  - abz(2,3)*difz(nz3)/difz(nz1)
      csx1(2,3)  = ksx(2,3) + asx(2,3)*difz(nz3)/difz(nz1)
      csx2(2,3)  =  - asx(2,3)*difz(nz3)/difz(nz1)
      csy1(2,3)  = ksy(2,3) + asy(2,3)*difz(nz3)/difz(nz1)
      csy2(2,3)  =  - asy(2,3)*difz(nz3)/difz(nz1)
      csz1(2,3)  = ksz(2,3) + asz(2,3)*difz(nz3)/difz(nz1)
      csz2(2,3)  =  - asz(2,3)*difz(nz3)/difz(nz1)
c
      do 10 j = 1,3
	if (perio(j)) then
	   lsym(1,j) = .false.
	   lsym(2,j) = .false.
	end if
  10  continue
c
      do 100 i = 1,2
       do 100 j = 1,3
	if ( perio(j) .or. lsym(i,j) ) then
	  crho1(i,j)  = 0.0
	  crho2(i,j)  = 0.0
	  cu1(i,j)    = 0.0
	  cu2(i,j)    = 0.0
	  cres1(i,j)  = 0.0
	  cres2(i,j)  = 0.0
	  cbx1(i,j)   = 0.0
	  cbx2(i,j)   = 0.0
	  cby1(i,j)   = 0.0
	  cby2(i,j)   = 0.0
	  cbz1(i,j)   = 0.0
	  cbz2(i,j)   = 0.0
	  csx1(i,j)   = 0.0
	  csx2(i,j)   = 0.0
	  csy1(i,j)   = 0.0
	  csy2(i,j)   = 0.0
	  csz1(i,j)   = 0.0
	  csz2(i,j)   = 0.0
	end if
  100 continue
      return
      end
