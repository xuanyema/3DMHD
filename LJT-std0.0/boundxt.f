	subroutine boundxt
c***********************************************************************
c
	include 'misflin'
c
	integer  ix,iz,iy
	real     rhomin,bymin,bzmin,sxmin,omega,umin,cot,sot,gaminv
c
c--------------------------------------
c   x = xmin
c
c--------------------------------------
        gaminv = 1./gamma
        omega=2*pi/20
        omega=omega*time
        cot=cos(omega)
        sot=sin(omega)       

        rhomin=rho0+delrho
         sxmin=rhomin*vx0        
         bymin=by0*cot+bz0*sot
         bzmin=bz0*cot-by0*sot

       do 110 iz = 2, nz1
       do 110 iy = 2, ny1
	 rho(1,iy,iz) = rhomin
	 sx(1,iy,iz)  = sxmin
	 sy(1,iy,iz)  = 0.0
	 sz(1,iy,iz)  = 0.0

	 by(1,iy,iz)  = bymin
	 bz(1,iy,iz)  = bzmin
c------div B=0.
	 bx(1,iy,iz)  = bx(3,iy,iz)+
     +            ( dify(iy)*(by(2,iy+1,iz)-by(2,iy-1,iz))
     +            + difz(iz)*(bz(2,iy,iz+1)-bz(2,iy,iz-1)) ) /difx(2)

c          u(1,iy,iz)= (bx(2,iy,iz)*bx(2,iy,iz)-bx(1,iy,iz)*bx(1,iy,iz))
c     +               +(by(2,iy,iz)*by(2,iy,iz)-by(1,iy,iz)*by(1,iy,iz))
c     +               +(bz(2,iy,iz)*bz(2,iy,iz)-bz(1,iy,iz)*bz(1,iy,iz))
c     +               + 2.*u(2,iy,iz)**gamma
c          if (u(1,iy,iz)<0.14) write(*,*) iy,iz, u(1,iy,iz)
c          u(1,iy,iz) = max(u(1,iy,iz),0.2)
          u(1,iy,iz) = pmsp

	  u(1,iy,iz) = (0.5*u(1,iy,iz) )**gaminv
	res(1,iy,iz) = res(3,iy,iz)
  110  continue

 


      end
