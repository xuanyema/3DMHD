	subroutine diag0
c***********************************************************************
	include 'misflin'
c
        integer  iz,ix,iy
c.....................................................
c   calculate energy density
c.....................................................
      do 40 iz = 1,nz
      do 40 iy = 1,ny
      do 40 ix = 1,nx
        help(ix,iy,iz)  = u(ix,iy,iz)**gamma
	heb(ix,iy,iz)   = 0.5*( bx(ix,iy,iz)*bx(ix,iy,iz)
     +                         + by(ix,iy,iz)*by(ix,iy,iz)
     +                          + bz(ix,iy,iz)*bz(ix,iy,iz) )
	hev(ix,iy,iz)   = 0.5*( sx(ix,iy,iz)*sx(ix,iy,iz)
     +                    + sy(ix,iy,iz)*sy(ix,iy,iz)
     +                     + sz(ix,iy,iz)*sz(ix,iy,iz) )/rho(ix,iy,iz) 
   40 continue
      return
      end
