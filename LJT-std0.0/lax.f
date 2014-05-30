	subroutine lax
c***********************************************************************
	include 'misflin'
c
c       Version 9d
c       This version is based on Version 8d.  It adds density
c       diagnostics to the magnetic field and pressure diagnostics
c       of Version 8d.
c       30 March 2004
c
	dimension term1(3),term2(3),jcb(3)
	integer  ix,iy,iz,ixanf
	integer  ix0,iy0,iz0
	real     eingam,chelp
	real     term1,term2,jcb
	real     prezok,hilpx,hilpy,hilpz
	real     rhilpx,rhilpy,rhilpz,rhelp
	logical  lndo
c .........................................................
c  integration of 3-d mhd equations
c .........................................................
      eingam = 1 - gamma
c--------------------------------------------------------------
c     Enhanced diagnostic variables
c     Determine whether or not enhanced diagnostics are generated.
      lndo = .false.
c
c     Define target coordinates
      ix0 = 62
      iy0 = 18
      iz0 = 39
c --------------------------------------------------------------
      do 20 iz = 1,nz
      do 20 iy = 1,ny
      ixanf = 1 + mod(iz+iy+1,2)
      do 20 ix = ixanf,nx,2	  
         hilf(ix,iy,iz) = 1./rho(ix,iy,iz)
   20 continue
      do 40 iz = 1,nz
      do 40 iy = 1,ny
      ixanf = 1 + mod(iz+iy+1,2)
      do 40 ix = ixanf,nx,2	  
         vx(ix,iy,iz) = sx(ix,iy,iz)*hilf(ix,iy,iz)
         vy(ix,iy,iz) = sy(ix,iy,iz)*hilf(ix,iy,iz)
         vz(ix,iy,iz) = sz(ix,iy,iz)*hilf(ix,iy,iz)
   40 continue
c
c************************************************
c    1. integration density
c    -------------------------
c************************************************
c
c     Density diagnostics
      if (lndo) then
c        Write density diagnostics before advancement
         write(6,27)
         write(6,27)
         write(6,*) 'Subroutine lax'
	 write(6,*) 'Density diagnostics'
         write(6,47) '    istep: ',istep
         write(6,48) '     time: ',time
         write(6,48) 'time step: ',dt
         write(6,27)
	 write(6,48) 'rho (before advancement): ',rho(ix0,iy0,iz0)
	 write(6,27)
	 rhilpx = difx(ix0)*(sx(ix0+1,iy0,iz0)-sx(ix0-1,iy0,iz0))
	 rhilpy = dify(iy0)*(sy(ix0,iy0+1,iz0)-sy(ix0,iy0-1,iz0))
	 rhilpz = difz(iz0)*(sz(ix0,iy0,iz0+1)-sz(ix0,iy0,iz0-1))
	 rhelp  = rhilpx + rhilpy + rhilpz
	 write(6,48) '   rhilpx: ',rhilpx
	 write(6,48) '   rhilpy: ',rhilpy
	 write(6,48) '   rhilpz: ',rhilpz
      endif
c
      do 100 iz = 2,nz-1
      do 100 iy = 2,ny1
      ixanf = 2 + mod(iz+iy+1,2)
      do 100 ix = ixanf,nx1,2
	  rho(ix,iy,iz) = rho(ix,iy,iz) - dt * (
     +                    difx(ix)*(sx(ix+1,iy,iz)-sx(ix-1,iy,iz))
     +                  + dify(iy)*(sy(ix,iy+1,iz)-sy(ix,iy-1,iz))
     +                  + difz(iz)*(sz(ix,iy,iz+1)-sz(ix,iy,iz-1)) )
  100 continue
c
      if (lndo) then
	 write(6,48) '    rhelp: ',rhelp
	 write(6,27)
	 write(6,48) ' rho (after advancement): ',rho(ix0,iy0,iz0)
	 write(6,27)
	 write(6,27)
      endif

c
c************************************************
c    2. momentum equations
c    -----------------------
c************************************************
c
      do 220 iz = 1,nz
      do 220 iy = 1,ny
      ixanf = 1 + mod(iz+iy+1,2)
      do 220 ix = ixanf,nx,2
	  help(ix,iy,iz) = u(ix,iy,iz)**gamma + 0.5*(
     +                         + bx(ix,iy,iz)*bx(ix,iy,iz)
     +                         + by(ix,iy,iz)*by(ix,iy,iz)
     +                         + bz(ix,iy,iz)*bz(ix,iy,iz) )
  220 continue
c
      do 240 iz = 1,nz
      do 240 iy = 1,ny
      ixanf = 1 + mod(iz+iy+1,2)
      do 240 ix = ixanf,nx,2
	  hilfx(ix,iy,iz) = help(ix,iy,iz)
     +                    + sx(ix,iy,iz)*vx(ix,iy,iz)
     +                    - bx(ix,iy,iz)*bx(ix,iy,iz)
	  hilfy(ix,iy,iz) = help(ix,iy,iz)
     +                    + sy(ix,iy,iz)*vy(ix,iy,iz)
     +                    - by(ix,iy,iz)*by(ix,iy,iz)
	  hilfz(ix,iy,iz) = help(ix,iy,iz)
     +                    + sz(ix,iy,iz)*vz(ix,iy,iz)
     +                    - bz(ix,iy,iz)*bz(ix,iy,iz)
	  feldx(ix,iy,iz) = sy(ix,iy,iz)*vz(ix,iy,iz)
     +                     - by(ix,iy,iz)*bz(ix,iy,iz)
	  feldy(ix,iy,iz) = sz(ix,iy,iz)*vx(ix,iy,iz)
     +                     - bz(ix,iy,iz)*bx(ix,iy,iz)
	  feldz(ix,iy,iz) = sx(ix,iy,iz)*vy(ix,iy,iz)
     +                     - bx(ix,iy,iz)*by(ix,iy,iz)
  240 continue
c
      do 260 iz = 2,nz-1
      do 260 iy = 2,ny1
      ixanf = 2 + mod(iz+iy+1,2)
      do 260 ix = ixanf,nx1,2
	  sx(ix,iy,iz) = sx(ix,iy,iz) - dt*(
     +              difx(ix)*(hilfx(ix+1,iy,iz)-hilfx(ix-1,iy,iz))
     +            + dify(iy)*(feldz(ix,iy+1,iz)-feldz(ix,iy-1,iz))
     +            + difz(iz)*(feldy(ix,iy,iz+1)-feldy(ix,iy,iz-1)) 
     +                 + muprof(ix,iy,iz)*sx(ix,iy,iz) )
	  sy(ix,iy,iz) = sy(ix,iy,iz) - dt*(
     +              dify(iy)*(hilfy(ix,iy+1,iz)-hilfy(ix,iy-1,iz))
     +            + difz(iz)*(feldx(ix,iy,iz+1)-feldx(ix,iy,iz-1))
     +            + difx(ix)*(feldz(ix+1,iy,iz)-feldz(ix-1,iy,iz)) 
     +                 + muprof(ix,iy,iz)*sy(ix,iy,iz) )
	  sz(ix,iy,iz) = sz(ix,iy,iz) - dt*(
     +              difz(iz)*(hilfz(ix,iy,iz+1)-hilfz(ix,iy,iz-1))
     +            + difx(ix)*(feldy(ix+1,iy,iz)-feldy(ix-1,iy,iz))
     +            + dify(iy)*(feldx(ix,iy+1,iz)-feldx(ix,iy-1,iz)) 
     +                 + muprof(ix,iy,iz)*sz(ix,iy,iz) )
  260 continue
c************************************************
c    3.a resistivity
c    ----------------
c************************************************
      do 499 iz = 1,nz
      do 499 iy = 1,ny
      do 499 ix = 1,nx
	  feldx(ix,iy,iz) = 0.
	  feldy(ix,iy,iz) = 0.
	  feldz(ix,iy,iz) = 0.
	  hilf(ix,iy,iz) = 0.
  499 continue

      do 501 iz = 2,nz-1
      do 501 iy = 2,ny1
      ixanf = 2 + mod(iz+iy+1,2)
      do 501 ix = ixanf,nx1,2
	  feldx(ix,iy,iz) = dify(iy)*(bz(ix,iy+1,iz)-bz(ix,iy-1,iz))
     +                   - difz(iz)*(by(ix,iy,iz+1)-by(ix,iy,iz-1))
	  feldy(ix,iy,iz) = difz(iz)*(bx(ix,iy,iz+1)-bx(ix,iy,iz-1))
     +                   - difx(ix)*(bz(ix+1,iy,iz)-bz(ix-1,iy,iz))
	  feldz(ix,iy,iz) = difx(ix)*(by(ix+1,iy,iz)-by(ix-1,iy,iz))
     +                   - dify(iy)*(bx(ix,iy+1,iz)-bx(ix,iy-1,iz))
	  hilf(ix,iy,iz) = feldx(ix,iy,iz)*feldx(ix,iy,iz)
     +                  + feldy(ix,iy,iz)*feldy(ix,iy,iz)
     +                  + feldz(ix,iy,iz)*feldz(ix,iy,iz)
  501   continue
      if (etasw.ne.0) call resist 
c************************************************
c    3. magnetic field
c    ------------------
c************************************************
c
c      hilf = v x b
c
      do 520 iz = 1,nz
      do 520 iy = 1,ny
      ixanf = 1 + mod(iz+iy+1,2)
      do 520 ix = ixanf,nx,2
	  hilfx(ix,iy,iz) = vy(ix,iy,iz)*bz(ix,iy,iz)
     +                      - vz(ix,iy,iz)*by(ix,iy,iz) 
	  hilfy(ix,iy,iz) = vz(ix,iy,iz)*bx(ix,iy,iz)
     +                      - vx(ix,iy,iz)*bz(ix,iy,iz) 
	  hilfz(ix,iy,iz) = vx(ix,iy,iz)*by(ix,iy,iz)
     +                      - vy(ix,iy,iz)*bx(ix,iy,iz)
  520 continue
c
c      vorsicht: dx,dy,dz -> res
c             ist abhaengig vom gewaehlten resitivitaetsmodell
c             (insbesondere fuer res =res(j) aendern!)
c
c     Enhanced diagnostics
      if (lndo) then
c        Write magnetic field diagnostics before advancement
         write(6,27)
         write(6,27)
         write(6,*) 'Subroutine lax'
	 write(6,*) 'Magnetic field diagnostics'
         write(6,47) '    istep: ',istep
         write(6,48) '     time: ',time
         write(6,48) 'time step: ',dt
         write(6,27)
         write(6,*) 'Central point'
         write(6,50) ix0,iy0,iz0
         write(6,55) x(ix0),y(iy0),z(iz0)
         write(6,27)
         write(6,48) '  difx: ',difx(ix0)
         write(6,48) ' ddifx:',ddifx(ix0)
         write(6,48) 'ddifmx:',ddifmx(ix0)
         write(6,48) 'ddifpx:',ddifpx(ix0)
         write(6,27)
         write(6,48) '  dify: ',dify(iy0)
         write(6,48) ' ddify:',ddify(iy0)
         write(6,48) 'ddifmy:',ddifmy(iy0)
         write(6,48) 'ddifpy:',ddifpy(iy0)
         write(6,27)
         write(6,48) '  difz: ',difz(iz0)
         write(6,48) ' ddifz:',ddifz(iz0)
         write(6,48) 'ddifmz:',ddifmz(iz0)
         write(6,48) 'ddifpz:',ddifpz(iz0)
         write(6,27)
         write(6,*) 'Magnetic field components'
         write(6,48) '  bx: ',bx(ix0,iy0,iz0)
         write(6,48) '  by: ',by(ix0,iy0,iz0)
         write(6,48) '  bz: ',bz(ix0,iy0,iz0)
         write(6,27)
         write(6,*) 'Current density'
         write(6,48) 'feldx:',feldx(ix0,iy0,iz0)
         write(6,48) 'feldy:',feldy(ix0,iy0,iz0)
         write(6,48) 'feldz:',feldz(ix0,iy0,iz0)
         write(6,27)
         write(6,*) 'Auxillary quantities'
         write(6,48) ' res: ',res(ix0,iy0,iz0)
         write(6,48) '   u: ',u(ix0,iy0,iz0)
         write(6,48) ' rho: ',rho(ix0,iy0,iz0)
         write(6,48) '  vx: ',vx(ix0,iy0,iz0)
         write(6,48) '  vy: ',vy(ix0,iy0,iz0)
         write(6,48) '  vz: ',vz(ix0,iy0,iz0)
         write(6,48) 'help: ',help(ix0,iy0,iz0)
         write(6,27)
         write(6,27)
         write(6,*) 'First neighbor in x'
         write(6,50) ix0-1,iy0,iz0
         write(6,55) x(ix0-1),y(iy0),z(iz0)
         write(6,48) '    u: ',u(ix0-1,iy0,iz0)
         write(6,48) '  rho: ',rho(ix0-1,iy0,iz0)
         write(6,48) '   vx: ',vx(ix0-1,iy0,iz0)
         write(6,48) '   vy: ',vy(ix0-1,iy0,iz0)
         write(6,48) '   vz: ',vz(ix0-1,iy0,iz0)
         write(6,48) '   bx: ',bx(ix0-1,iy0,iz0)
         write(6,48) '   by: ',by(ix0-1,iy0,iz0)
         write(6,48) '   bz: ',bz(ix0-1,iy0,iz0)
         write(6,48) 'hilfx: ',hilfx(ix0-1,iy0,iz0)
         write(6,48) 'hilfy: ',hilfy(ix0-1,iy0,iz0)
         write(6,48) 'hilfz: ',hilfz(ix0-1,iy0,iz0)
         write(6,48) '  res: ',res(ix0-1,iy0,iz0)
         write(6,27)
         write(6,27)
         write(6,*) 'Second neighbor in x'
         write(6,50) ix0+1,iy0,iz0
         write(6,55) x(ix0+1),y(iy0),z(iz0)
         write(6,48) '    u: ',u(ix0+1,iy0,iz0)
         write(6,48) '  rho: ',rho(ix0+1,iy0,iz0)
         write(6,48) '   vx: ',vx(ix0+1,iy0,iz0)
         write(6,48) '   vy: ',vy(ix0+1,iy0,iz0)
         write(6,48) '   vz: ',vz(ix0+1,iy0,iz0)
         write(6,48) '   bx: ',bx(ix0+1,iy0,iz0)
         write(6,48) '   by: ',by(ix0+1,iy0,iz0)
         write(6,48) '   bz: ',bz(ix0+1,iy0,iz0)
         write(6,48) 'hilfx: ',hilfx(ix0+1,iy0,iz0)
         write(6,48) 'hilfy: ',hilfy(ix0+1,iy0,iz0)
         write(6,48) 'hilfz: ',hilfz(ix0+1,iy0,iz0)
         write(6,48) '  res: ',res(ix0+1,iy0,iz0)
         write(6,27)
         write(6,27)
         write(6,*) 'First neighbor in y'
         write(6,50) ix0,iy0-1,iz0
         write(6,55) x(ix0),y(iy0-1),z(iz0)
         write(6,48) '    u: ',u(ix0,iy0-1,iz0)
         write(6,48) '  rho: ',rho(ix0,iy0-1,iz0)
         write(6,48) '   vx: ',vx(ix0,iy0-1,iz0)
         write(6,48) '   vy: ',vy(ix0,iy0-1,iz0)
         write(6,48) '   vz: ',vz(ix0,iy0-1,iz0)
         write(6,48) '   bx: ',bx(ix0,iy0-1,iz0)
         write(6,48) '   by: ',by(ix0,iy0-1,iz0)
         write(6,48) '   bz: ',bz(ix0,iy0-1,iz0)
         write(6,48) 'hilfx: ',hilfx(ix0,iy0-1,iz0)
         write(6,48) 'hilfy: ',hilfy(ix0,iy0-1,iz0)
         write(6,48) 'hilfz: ',hilfz(ix0,iy0-1,iz0)
         write(6,48) '  res: ',res(ix0,iy0-1,iz0)
         write(6,27)
         write(6,27)
         write(6,*) 'Second neighbor in y'
         write(6,50) ix0,iy0+1,iz0
         write(6,55) x(ix0),y(iy0+1),z(iz0)
         write(6,48) '    u: ',u(ix0,iy0+1,iz0)
         write(6,48) '  rho: ',rho(ix0,iy0+1,iz0)
         write(6,48) '   vx: ',vx(ix0,iy0+1,iz0)
         write(6,48) '   vy: ',vy(ix0,iy0+1,iz0)
         write(6,48) '   vz: ',vz(ix0,iy0+1,iz0)
         write(6,48) '   bx: ',bx(ix0,iy0+1,iz0)
         write(6,48) '   by: ',by(ix0,iy0+1,iz0)
         write(6,48) '   bz: ',bz(ix0,iy0+1,iz0)
         write(6,48) 'hilfx: ',hilfx(ix0,iy0+1,iz0)
         write(6,48) 'hilfy: ',hilfy(ix0,iy0+1,iz0)
         write(6,48) 'hilfz: ',hilfz(ix0,iy0+1,iz0)
         write(6,48) '  res: ',res(ix0,iy0+1,iz0)
         write(6,27)
         write(6,27)
         write(6,*) 'First neighbor in z'
         write(6,50) ix0,iy0,iz0-1
         write(6,55) x(ix0),y(iy0),z(iz0-1)
         write(6,48) '    u: ',u(ix0,iy0,iz0-1)
         write(6,48) '  rho: ',rho(ix0,iy0,iz0-1)
         write(6,48) '   vx: ',vx(ix0,iy0,iz0-1)
         write(6,48) '   vy: ',vy(ix0,iy0,iz0-1)
         write(6,48) '   vz: ',vz(ix0,iy0,iz0-1)
         write(6,48) '   bx: ',bx(ix0,iy0,iz0-1)
         write(6,48) '   by: ',by(ix0,iy0,iz0-1)
         write(6,48) '   bz: ',bz(ix0,iy0,iz0-1)
         write(6,48) 'hilfx: ',hilfx(ix0,iy0,iz0-1)
         write(6,48) 'hilfy: ',hilfy(ix0,iy0,iz0-1)
         write(6,48) 'hilfz: ',hilfz(ix0,iy0,iz0-1)
         write(6,48) '  res: ',res(ix0,iy0,iz0-1)
         write(6,27)
         write(6,27)
         write(6,*) 'Second neighbor in z'
         write(6,50) ix0,iy0,iz0+1
         write(6,55) x(ix0),y(iy0),z(iz0+1)
         write(6,48) '    u: ',u(ix0,iy0,iz0+1)
         write(6,48) '  rho: ',rho(ix0,iy0,iz0+1)
         write(6,48) '   vx: ',vx(ix0,iy0,iz0+1)
         write(6,48) '   vy: ',vy(ix0,iy0,iz0+1)
         write(6,48) '   vz: ',vz(ix0,iy0,iz0+1)
         write(6,48) '   bx: ',bx(ix0,iy0,iz0+1)
         write(6,48) '   by: ',by(ix0,iy0,iz0+1)
         write(6,48) '   bz: ',bz(ix0,iy0,iz0+1)
         write(6,48) 'hilfx: ',hilfx(ix0,iy0,iz0+1)
         write(6,48) 'hilfy: ',hilfy(ix0,iy0,iz0+1)
         write(6,48) 'hilfz: ',hilfz(ix0,iy0,iz0+1)
         write(6,48) '  res: ',res(ix0,iy0,iz0+1)
         write(6,27)
         write(6,27)
         write(6,27)
c
c        Define variables corresponding to different terms in the 
c        equations advancing the magnetic field.
c        term1: the curl of v cross b
c        term2: j cross grad eta
c        Note that term3 (eta Laplacian b) does not appear in this
c        subroutine (although it appears in leap).
c        Array indices 1, 2, and 3 correspond to x, y, and z components
c        of b to which the term belongs
c
         term1(1) = dify(iy0)*( 
     +                   hilfz(ix0,iy0+1,iz0)-hilfz(ix0,iy0-1,iz0) )
     +            - difz(iz0)*( 
     +                   hilfy(ix0,iy0,iz0+1)-hilfy(ix0,iy0,iz0-1) )
         term2(1) = 
     +        - feldz(ix0,iy0,iz0)*dify(iy0)*(
     +                        res(ix0,iy0+1,iz0)-res(ix0,iy0-1,iz0) )     
     +        + feldy(ix0,iy0,iz0)*difz(iz0)*(
     +                        res(ix0,iy0,iz0+1)-res(ix0,iy0,iz0-1) )
c
         term1(2) = 
     +             difz(iz0)*( 
     +                  hilfx(ix0,iy0,iz0+1)-hilfx(ix0,iy0,iz0-1) ) 
     +           - difx(ix0)*( 
     +                  hilfz(ix0+1,iy0,iz0)-hilfz(ix0-1,iy0,iz0) )
         term2(2) =
     +        - feldx(ix0,iy0,iz0)*difz(iz0)*(
     +                        res(ix0,iy0,iz0+1)-res(ix0,iy0,iz0-1))
     +       + feldz(ix0,iy0,iz0)*difx(ix0)*(
     +                        res(ix0+1,iy0,iz0)-res(ix0-1,iy0,iz0))
c
         term1(3) =
     +             difx(ix0)*( 
     +                  hilfy(ix0+1,iy0,iz0)-hilfy(ix0-1,iy0,iz0) )
     +           - dify(iy0)*( 
     +                  hilfx(ix0,iy0+1,iz0)-hilfx(ix0,iy0-1,iz0) )
         term2(3) =
     +        - feldy(ix0,iy0,iz0)*difx(ix0)*(
     +                        res(ix0+1,iy0,iz0)-res(ix0-1,iy0,iz0))
     +        + feldx(ix0,iy0,iz0)*dify(iy0)*(
     +                        res(ix0,iy0+1,iz0)-res(ix0,iy0-1,iz0))
c
c        Components of j cross b
         jcb(1) = feldy(ix0,iy0,iz0)*bz(ix0,iy0,iz0) 
     +          - feldz(ix0,iy0,iz0)*by(ix0,iy0,iz0)
         jcb(2) = feldz(ix0,iy0,iz0)*bx(ix0,iy0,iz0)
     +          - feldx(ix0,iy0,iz0)*bz(ix0,iy0,iz0)
         jcb(3) = feldx(ix0,iy0,iz0)*by(ix0,iy0,iz0)
     +          - feldy(ix0,iy0,iz0)*bx(ix0,iy0,iz0)
c
c
c        Report the values of the various terms.
         write(6,*) 'Advancement terms'
         write(6,*) '-----------------'
         write(6,27)
         write(6,27)
         write(6,*) 'bx'
         write(6,*) '--'
         write(6,48) 'curl ( v X b ):',term1(1)
         write(6,48) '- grad ( eta ) X j:',term2(1)
         write(6,27)
         write(6,27)
         write(6,*) 'by'
         write(6,*) '--'
         write(6,48) 'curl ( v X b ):',term1(2)
         write(6,48) '- grad ( eta ) X j:',term2(2)
         write(6,27)
         write(6,27)
         write(6,*) 'bz'
         write(6,*) '--'
         write(6,48) 'curl ( v X b ):',term1(3)
         write(6,48) '- grad ( eta ) X j:',term2(3)
         write(6,27)
         write(6,27)
         write(6,27)
         write(6,*) 'j X b components'
         write(6,*) '----------------'
         write(6,48) '(j X b)_x:',jcb(1)
         write(6,48) '(j X b)_y:',jcb(2)
         write(6,48) '(j X b)_z:',jcb(3)
         write(6,27)
         write(6,27)
         write(6,27)
      endif
c
c
c
c     The actual advancement of the magnetic field
      do 560 iz = 2,nz-1
      do 560 iy = 2,ny1
      ixanf = 2 + mod(iz+iy+1,2)
      do 560 ix = ixanf,nx1,2
	  bx(ix,iy,iz) = bx(ix,iy,iz) + dt*(
     +            dify(iy)*( hilfz(ix,iy+1,iz)-hilfz(ix,iy-1,iz) )
     +          - difz(iz)*( hilfy(ix,iy,iz+1)-hilfy(ix,iy,iz-1) )
     +    - feldz(ix,iy,iz)*dify(iy)*(res(ix,iy+1,iz)-res(ix,iy-1,iz)) 
     +    + feldy(ix,iy,iz)*difz(iz)*(res(ix,iy,iz+1)-res(ix,iy,iz-1)) )
	  by(ix,iy,iz) = by(ix,iy,iz) + dt*(
     +            difz(iz)*( hilfx(ix,iy,iz+1)-hilfx(ix,iy,iz-1) )
     +          - difx(ix)*( hilfz(ix+1,iy,iz)-hilfz(ix-1,iy,iz) )
     +    - feldx(ix,iy,iz)*difz(iz)*(res(ix,iy,iz+1)-res(ix,iy,iz-1))
     +    + feldz(ix,iy,iz)*difx(ix)*(res(ix+1,iy,iz)-res(ix-1,iy,iz)) )
	  bz(ix,iy,iz) = bz(ix,iy,iz) + dt*(
     +            difx(ix)*( hilfy(ix+1,iy,iz)-hilfy(ix-1,iy,iz) )
     +          - dify(iy)*( hilfx(ix,iy+1,iz)-hilfx(ix,iy-1,iz) ) 
     +    - feldy(ix,iy,iz)*difx(ix)*(res(ix+1,iy,iz)-res(ix-1,iy,iz))
     +    + feldx(ix,iy,iz)*dify(iy)*(res(ix,iy+1,iz)-res(ix,iy-1,iz)) )
  560 continue
c
      if (lndo) then
	 write(6,*) 'After advancement'
	 write(6,*) '-----------------'
	 write(6,27)
	 write(6,48) 'bx: ',bx(ix0,iy0,iz0)
	 write(6,48) 'by: ',by(ix0,iy0,iz0)
	 write(6,48) 'bz: ',bz(ix0,iy0,iz0)
	 write(6,27)
	 write(6,27)
	 write(6,27)
      endif
c
c************************************************
c    4. pressure
c    --------------
c************************************************
c
c      p = 2 * u**gamma
c
c
      do 820 iz = 1,nz
      do 820 iy = 1,ny
      ixanf = 1 + mod(iz+iy+1,2)
      do 820 ix = ixanf,nx,2
	  hilfx(ix,iy,iz) = u(ix,iy,iz)*vx(ix,iy,iz)
	  hilfy(ix,iy,iz) = u(ix,iy,iz)*vy(ix,iy,iz)
	  hilfz(ix,iy,iz) = u(ix,iy,iz)*vz(ix,iy,iz)
  820 continue
c
c
c     Additional pressure diagnostics
      if (lndo) then
	 write(6,27)
	 write(6,27)
	 write(6,*) 'Subroutine lax'
	 write(6,*) 'Pressure diagnostics'
	 write(6,47) 'istep: ',istep
	 write(6,48) ' time: ',time
	 write(6,48) 'time step: ',dt
	 write(6,27)
	 write(6,*) 'Central point'
	 write(6,50) ix0,iy0,iz0
	 write(6,55) x(ix0),y(iy0),z(iz0)
	 write(6,48) 'difx: ',difx(ix0)
	 write(6,48) 'dify: ',dify(iy0)
	 write(6,48) 'difz: ',difz(iz0)
	 write(6,48) ' res: ',res(ix0,iy0,iz0)
 	 write(6,48) ' rho: ',rho(ix0,iy0,iz0)
	 write(6,48) '  vx: ',vx(ix0,iy0,iz0)
	 write(6,48) '  vy: ',vy(ix0,iy0,iz0)
	 write(6,48) '  vz: ',vz(ix0,iy0,iz0)
c	 write(6,48) 'hilf: ',hilf(ix0,iy0,iz0)
	 write(6,48) 'u (before advancement): ',u(ix0,iy0,iz0)
	 write(6,27)
	 write(6,27)
	 write(6,*) 'First neighbor in x'
	 write(6,50) ix0-1,iy0,iz0
	 write(6,55) x(ix0-1),y(iy0),z(iz0)
	 write(6,48) '    u: ',u(ix0-1,iy0,iz0)
	 write(6,48) '   vx: ',vx(ix0-1,iy0,iz0)
	 write(6,48) '   vy: ',vy(ix0-1,iy0,iz0)
	 write(6,48) '   vz: ',vz(ix0-1,iy0,iz0)
	 write(6,48) 'hilfx: ',hilfx(ix0-1,iy0,iz0)
	 write(6,27)
	 write(6,*) 'Second neighbor in x'
	 write(6,50) ix0+1,iy0,iz0
	 write(6,55) x(ix0+1),y(iy0),z(iz0)
	 write(6,48) '    u: ',u(ix0+1,iy0,iz0)
	 write(6,48) '   vx: ',vx(ix0+1,iy0,iz0)
	 write(6,48) '   vy: ',vy(ix0+1,iy0,iz0)
	 write(6,48) '   vz: ',vz(ix0+1,iy0,iz0)
	 write(6,48) 'hilfx: ',hilfx(ix0+1,iy0,iz0)
	 write(6,27)
	 write(6,27)
	 write(6,*) 'First neighbor in y'
	 write(6,50) ix0,iy0-1,iz0
	 write(6,55) x(ix0),y(iy0-1),z(iz0)
	 write(6,48) '    u: ',u(ix0,iy0-1,iz0)
	 write(6,48) '   vx: ',vx(ix0,iy0-1,iz0)
	 write(6,48) '   vy: ',vy(ix0,iy0-1,iz0)
	 write(6,48) '   vz: ',vz(ix0,iy0-1,iz0)
	 write(6,48) 'hilfy: ',hilfy(ix0,iy0-1,iz0)
	 write(6,27)
	 write(6,27)
	 write(6,*) 'Second neighbor in y'
	 write(6,50) ix0,iy0+1,iz0
	 write(6,55) x(ix0),y(iy0+1),z(iz0)
	 write(6,48) '    u: ',u(ix0,iy0+1,iz0)
	 write(6,48) '   vx: ',vx(ix0,iy0+1,iz0)
	 write(6,48) '   vy: ',vy(ix0,iy0+1,iz0)
	 write(6,48) '   vz: ',vz(ix0,iy0+1,iz0)
	 write(6,48) 'hilfy: ',hilfy(ix0,iy0+1,iz0)
	 write(6,27)
	 write(6,27)
	 write(6,*) 'First neighbor in z'
	 write(6,50) ix0,iy0,iz0-1
	 write(6,55) x(ix0),y(iy0),z(iz0-1)
	 write(6,48) ' u: ',u(ix0,iy0,iz0-1)
	 write(6,48) 'vx: ',vx(ix0,iy0,iz0-1)
	 write(6,48) 'vy: ',vy(ix0,iy0,iz0-1)
	 write(6,48) 'vz: ',vz(ix0,iy0,iz0-1)
	 write(6,48) 'hilfz: ',hilfz(ix0,iy0,iz0-1)
	 write(6,27)
	 write(6,27)
	 write(6,*) 'Second neighbor in z'
	 write(6,50) ix0,iy0,iz0+1
	 write(6,55) x(ix0),y(iy0),z(iz0+1)
	 write(6,48) ' u: ',u(ix0,iy0,iz0+1)
	 write(6,48) 'vx: ',vx(ix0,iy0,iz0+1)
	 write(6,48) 'vy: ',vy(ix0,iy0,iz0+1)
	 write(6,48) 'vz: ',vz(ix0,iy0,iz0+1)
	 write(6,48) 'hilfz: ',hilfz(ix0,iy0,iz0+1)
	 write(6,27)
	 write(6,27)
	 write(6,27)
	 hilpx = difx(ix0)*(hilfx(ix0+1,iy0,iz0)-hilfx(ix0-1,iy0,iz0))
	 hilpy = dify(iy0)*(hilfy(ix0,iy0+1,iz0)-hilfy(ix0,iy0-1,iz0))
	 hilpz = difz(iz0)*(hilfz(ix0,iy0,iz0+1)-hilfz(ix0,iy0,iz0-1))
	 write(6,*) 'Derived quantities'
c	 write(6,48) ' help: ',help(ix0,iy0,iz0)
	 write(6,48) 'hilpx: ',hilpx
	 write(6,48) 'hilpy: ',hilpy
	 write(6,48) 'hilpz: ',hilpz
	 write(6,27)
	 prezok = hilpx + hilpy + hilpz
	 write(6,48) ' Sum of hilp terms: ',prezok
	 write(6,48) '  Advancement term: ',dt*prezok
	 write(6,48) 'Estimate for new u: ',u(ix0,iy0,iz0) - dt*prezok
	 write(6,27)
      endif
c
c
c

      chelp = eingam/gamma/3.0**eingam
      do 840 iz = 2,nz-1
      do 840 iy = 2,ny1
      ixanf = 2 + mod(iz+iy+1,2)
      do 840 ix = ixanf,nx1,2
	  u(ix,iy,iz) = u(ix,iy,iz) - dt*(
     +            difx(ix)*(hilfx(ix+1,iy,iz)-hilfx(ix-1,iy,iz))
     +          + dify(iy)*(hilfy(ix,iy+1,iz)-hilfy(ix,iy-1,iz))
     +          + difz(iz)*(hilfz(ix,iy,iz+1)-hilfz(ix,iy,iz-1))
     +          + chelp*res(ix,iy,iz)*hilf(ix,iy,iz)*(
     +            ( meanpx(ix)*u(ix+1,iy,iz)+meanmx(ix)*u(ix-1,iy,iz) )
     +           +( meanpy(iy)*u(ix,iy+1,iz)+meanmy(iy)*u(ix,iy-1,iz) )
     +           +( meanpz(iz)*u(ix,iy,iz+1)+meanmz(iz)*u(ix,iy,iz-1) )
     +                           )**eingam )
  840 continue
c
      if ((mod((istep+1), nrelax) .eq. 0) .and. relax 
     +        .and.   (istep+1 .le. nmrelax) ) then
	 do 937 iz = 1,nz
	 do 937 iy = 1,ny
	 do 937 ix = 1,nx
	   sx(ix,iy,iz) = 0.0
	   sy(ix,iy,iz) = 0.0
	   sz(ix,iy,iz) = 0.0
 937	continue
      endif
c
c----------------------------------------------------------------------
c

  27  format((/))
  45  format(1x,f7.2, 4f10.4)
  47  format( 1x,A,i5)
  48  format( 1x,A,f14.5)
  50  format(1x,'Indices: (',i3,',',i3,',',i3,')')
  55  format(1x,'Location: (',f8.3,',',f8.3,',',f8.3,')')


      return
      end
