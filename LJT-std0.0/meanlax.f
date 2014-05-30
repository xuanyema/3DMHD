	subroutine meanlax
c***********************************************************************
	include 'misflin'
c
	integer  ix,iy,iz,ixanf
c .................................................................
c  averages for the lax integration
c .................................................................

      do 11 iz = 2,nz1
       do 11 iy = 2,ny1
	ixanf = 2 + mod(iz+iy+1,2)
	do 11 ix = ixanf,nx1,2
	  rho(ix,iy,iz)  = ( (meanpx(ix)*rho(ix+1,iy,iz)
     +                       + meanmx(ix)*rho(ix-1,iy,iz) )
     +                  + (meanpy(iy)*rho(ix,iy+1,iz)
     +                       + meanmy(iy)*rho(ix,iy-1,iz) )
     +                  + (meanpz(iz)*rho(ix,iy,iz+1)
     +                       + meanmz(iz)*rho(ix,iy,iz-1) ) )/3.
	  u(ix,iy,iz)    = ( (meanpx(ix)*u(ix+1,iy,iz)
     +                       + meanmx(ix)*u(ix-1,iy,iz) )
     +                  + (meanpy(iy)*u(ix,iy+1,iz)
     +                       + meanmy(iy)*u(ix,iy-1,iz) )
     +                  + (meanpz(iz)*u(ix,iy,iz+1)
     +                       + meanmz(iz)*u(ix,iy,iz-1) ) )/3.
	  bx(ix,iy,iz)  = ( (meanpx(ix)*bx(ix+1,iy,iz)
     +                       + meanmx(ix)*bx(ix-1,iy,iz) )
     +                  + (meanpy(iy)*bx(ix,iy+1,iz)
     +                       + meanmy(iy)*bx(ix,iy-1,iz) )
     +                  + (meanpz(iz)*bx(ix,iy,iz+1)
     +                       + meanmz(iz)*bx(ix,iy,iz-1) ) )/3.
	  by(ix,iy,iz)  = ( (meanpx(ix)*by(ix+1,iy,iz)
     +                       + meanmx(ix)*by(ix-1,iy,iz) )
     +                  + (meanpy(iy)*by(ix,iy+1,iz)
     +                       + meanmy(iy)*by(ix,iy-1,iz) )
     +                  + (meanpz(iz)*by(ix,iy,iz+1)
     +                       + meanmz(iz)*by(ix,iy,iz-1) ) )/3.
	  bz(ix,iy,iz)  = ( (meanpx(ix)*bz(ix+1,iy,iz)
     +                       + meanmx(ix)*bz(ix-1,iy,iz) )
     +                  + (meanpy(iy)*bz(ix,iy+1,iz)
     +                       + meanmy(iy)*bz(ix,iy-1,iz) )
     +                  + (meanpz(iz)*bz(ix,iy,iz+1)
     +                       + meanmz(iz)*bz(ix,iy,iz-1) ) )/3.
	  sx(ix,iy,iz)  = ( (meanpx(ix)*sx(ix+1,iy,iz)
     +                       + meanmx(ix)*sx(ix-1,iy,iz) )
     +                  + (meanpy(iy)*sx(ix,iy+1,iz)
     +                       + meanmy(iy)*sx(ix,iy-1,iz) )
     +                  + (meanpz(iz)*sx(ix,iy,iz+1)
     +                       + meanmz(iz)*sx(ix,iy,iz-1) ) )/3.
	  sy(ix,iy,iz)  = ( (meanpx(ix)*sy(ix+1,iy,iz)
     +                       + meanmx(ix)*sy(ix-1,iy,iz) )
     +                  + (meanpy(iy)*sy(ix,iy+1,iz)
     +                       + meanmy(iy)*sy(ix,iy-1,iz) )
     +                  + (meanpz(iz)*sy(ix,iy,iz+1)
     +                       + meanmz(iz)*sy(ix,iy,iz-1) ) )/3.
	  sz(ix,iy,iz)  = ( (meanpx(ix)*sz(ix+1,iy,iz)
     +                       + meanmx(ix)*sz(ix-1,iy,iz) )
     +                  + (meanpy(iy)*sz(ix,iy+1,iz)
     +                       + meanmy(iy)*sz(ix,iy-1,iz) )
     +                  + (meanpz(iz)*sz(ix,iy,iz+1)
     +                       + meanmz(iz)*sz(ix,iy,iz-1) ) )/3.
   11 continue
c
      return
      end
