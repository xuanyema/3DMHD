	subroutine smooth
c***********************************************************************
c
c       Version 3
c       This version simply formats the code that some of the logic is
c       easier to see upon casual inspection.
c       Again, Version 2 incorporates local smoothing, as instituted by
c       Antonius.  It should be used with Version 4 of termin and
c       Version 16 of misflin (and perhaps later versions of that
c       subroutine and that include file, as well).
c       24 February 2004
c
	include 'misflin'
c
	integer  ix,iy,iz
	real     ppp1(nz),ppp2(nz)
c ................................................................
	do 50 iz = 1,nz
	   ppp1(iz) =  0.25
           ppp2(iz) =  0.25
c	   ppp1(iz) =  (1.-0.75/cosh(4.*z(iz)/zmax))
c          ppp2(iz) =  0.25/cosh(4.*z(iz)/zmax)
   50   continue
c
        do 100 iz = 2,nz1
	do 100 iy = 2,ny1
	do 100 ix = 2,nx1
	   if (rhocor(ix,iy,iz)) then
	      help(ix,iy,iz) =  ppp1(iz)*rho(ix,iy,iz)
     +                  + ppp2(iz)*(meanpx(ix)*rho(ix+1,iy,iz)
     +                        + meanmx(ix)*rho(ix-1,iy,iz)
     +                        + meanpy(iy)*rho(ix,iy+1,iz)
     +                        + meanmy(iy)*rho(ix,iy-1,iz) 
     +                        + meanpz(iz)*rho(ix,iy,iz+1)
     +                        + meanmz(iz)*rho(ix,iy,iz-1) )
	   endif
  100   continue
c
        do 101 iz = 2,nz1
	do 101 iy = 2,ny1
	do 101 ix = 2,nx1
 	   if (rhocor(ix,iy,iz)) then
c             write (*,*) 'ix,iy,iz,rho: ',ix,iy,iz,rho(ix,iy,iz)
	      rho(ix,iy,iz) = help(ix,iy,iz)
c             write (*,*) 'rho: ',rho(ix,iy,iz)
	   endif
  101   continue
c
        do 200 iz = 2,nz1
	do 200 iy = 2,ny1
	do 200 ix = 2,nx1
	   if (ucor(ix,iy,iz)) then
	      help(ix,iy,iz) =  ppp1(iz)*u(ix,iy,iz)
     +              + ppp2(iz)*(meanpx(ix)*u(ix+1,iy,iz)
     +                        + meanmx(ix)*u(ix-1,iy,iz)
     +                        + meanpy(iy)*u(ix,iy+1,iz)
     +                        + meanmy(iy)*u(ix,iy-1,iz)
     +                        + meanpz(iz)*u(ix,iy,iz+1)
     +                        + meanmz(iz)*u(ix,iy,iz-1) )
	   endif
  200   continue
c
	
        do 201 iz = 2,nz1
        do 201 iy = 2,ny1
        do 201 ix = 2,nx1
	   if (ucor(ix,iy,iz)) then
	      u(ix,iy,iz) = help(ix,iy,iz)
	   endif
  201   continue
c
	
        do 300 iz = 2,nz1
        do 300 iy = 2,ny1
        do 300 ix = 2,nx1
           if (sxyzcor(ix,iy,iz)) then
	      feldx(ix,iy,iz) =  ppp1(iz)*sx(ix,iy,iz)
     +              + ppp2(iz)*(meanpx(ix)*sx(ix+1,iy,iz)
     +                        + meanmx(ix)*sx(ix-1,iy,iz)
     +                        + meanpy(iy)*sx(ix,iy+1,iz)
     +                        + meanmy(iy)*sx(ix,iy-1,iz)
     +                        + meanpz(iz)*sx(ix,iy,iz+1)
     +                        + meanmz(iz)*sx(ix,iy,iz-1) )
	      feldy(ix,iy,iz) =  ppp1(iz)*sy(ix,iy,iz)
     +              + ppp2(iz)*(meanpx(ix)*sy(ix+1,iy,iz)
     +                        + meanmx(ix)*sy(ix-1,iy,iz)
     +                        + meanpy(iy)*sy(ix,iy+1,iz)
     +                        + meanmy(iy)*sy(ix,iy-1,iz)
     +                        + meanpz(iz)*sy(ix,iy,iz+1)
     +                        + meanmz(iz)*sy(ix,iy,iz-1) )
	      feldz(ix,iy,iz) =  ppp1(iz)*sz(ix,iy,iz)
     +              + ppp2(iz)*(meanpx(ix)*sz(ix+1,iy,iz)
     +                        + meanmx(ix)*sz(ix-1,iy,iz)
     +                        + meanpy(iy)*sz(ix,iy+1,iz)
     +                        + meanmy(iy)*sz(ix,iy-1,iz)
     +                        + meanpz(iz)*sz(ix,iy,iz+1)
     +                        + meanmz(iz)*sz(ix,iy,iz-1) )
           endif
  300   continue
c
        do 301 iz = 2,nz1
        do 301 iy = 2,ny1
        do 301 ix = 2,nx1
	   if (sxyzcor(ix,iy,iz)) then
c             write (*,*) 'ix,iy,iz,sx,sy,sz: ',
c    +                   ix,iy,iz,
c    +                   sx(ix,iy,iz),sy(ix,iy,iz),sz(ix,iy,iz)
	      sx(ix,iy,iz) = feldx(ix,iy,iz)
	      sy(ix,iy,iz) = feldy(ix,iy,iz)
	      sz(ix,iy,iz) = feldz(ix,iy,iz)
c             write (*,*) 'ix,iy,iz,sx,sy,sz: ',
c    +                   ix,iy,iz,
c    +                   sx(ix,iy,iz),sy(ix,iy,iz),sz(ix,iy,iz)
           endif
  301   continue
c
c
        return
        end
