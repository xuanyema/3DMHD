	subroutine leap
c***********************************************************************
	include 'misflin'
c
c       Version 17d
c       Version 17d is based on Version 16d.  It adds density 
c       diagnostics to the magnetic field and pressure diagnostics
c       of Version 16d.
c       by Antonius Otto and Fred Hall IV
c       30 March 2004
c
	dimension term1(3),term2(3),term3(3),jcb(3)
	integer  ix,iy,iz,ixanf,ixx,iyy,izz
	integer  kx1,kx2,ky,kz
	integer  npdiag,rsw,newrelax,rsw1,newrelax1,
     +              rsw2,newrelax2,rsw3,newrelax3,rsw4,newrelax4,
     +              rsw5,newrelax5,rsw6,newrelax6,rsw7,newrelax7,
     +              rsw8,newrelax8
	integer  ix0,iy0,iz0
	real     ddt,eingam,dif0,chelp,pmin,pinf,gaminv,
     +           hilpx(nx,ny,nz),hilpy(nx,ny,nz),hilpz(nx,ny,nz)
	real     dfux,dfuy,dfuz
	real     mu,mu0,musw
	real     term1,term2,term3,dterm,jcb
	real     prezok,hlpstore
	real     rhilpx,rhilpy,rhilpz
	logical  lndo
c .............................................................
c  integration of 3-d mhdequations
c .............................................................
      ky=12
      kz=8
      npdiag = 200
      ddt = 2.*dt
      eingam = 1 - gamma
      gaminv = 1./gamma
c--------------------------------------------------------------
c     Set nrelax switch variables
      rsw = 400
      newrelax = 200
      rsw1 = 800
      newrelax1 = 400
      rsw2 = 1600
      newrelax2 = 800
      rsw3 = 8000
      newrelax3 = 80
      rsw4 = 8800
      newrelax4 = 800
      rsw5 = 14400
      newrelax5 = 80
      rsw6 = 15200
      newrelax6 = 1600
      rsw7 = 19200
      newrelax7 = 80
      rsw8 = 20000
      newrelax8 = 800
c-----------------------------------------------------------------
c     Enhanced diagnostic variables
c     Determine whether or not enhanced diagnostics are generated.
      lndo = .false.
c
c     Define target coordinates
      ix0 = 62
      iy0 = 18
      iz0 = 39
c -----------------------------------------------------------------
      if (zentr(1)) then
        dif0 = ddifx((nx+1)/2)
      else
        dif0 = ddifx(2)
      end if
      do 20 iz = 1,nzz
      do 20 iy = 1,ny
      ixanf = 1 + mod(iz+iy+igrid,2)
      do 20 ix = ixanf,nx,2	  
         hilf(ix,iy,iz) = 1./rho(ix,iy,iz)
   20 continue
      do 40 iz = 1,nzz
      do 40 iy = 1,ny
      ixanf = 1 + mod(iz+iy+igrid,2)
      do 40 ix = ixanf,nx,2	  
         vx(ix,iy,iz) = sx(ix,iy,iz)*hilf(ix,iy,iz)
         vy(ix,iy,iz) = sy(ix,iy,iz)*hilf(ix,iy,iz)
         vz(ix,iy,iz) = sz(ix,iy,iz)*hilf(ix,iy,iz)
   40 continue
      
c
c************************************************
c    1. density
c    ------------
c************************************************
c
c     Density diagnostics
      if (lndo) then
c        Write density diagnostics before advancement
         write(6,27)
         write(6,27)
         write(6,*) 'Subroutine leap'
	 write(6,*) 'Density diagnostics'
         write(6,47) '    istep: ',istep
         write(6,48) '     time: ',time
         write(6,48) 'time step: ',dt
         write(6,48) 'double dt: ',ddt
         write(6,27)
	 write(6,48) 'rho (before advancement): ',rho(ix0,iy0,iz0)
	 write(6,27)
	 rhilpx = difx(ix0)*(sx(ix0+1,iy0,iz0)-sx(ix0-1,iy0,iz0))
	 rhilpy = dify(iy0)*(sy(ix0,iy0+1,iz0)-sy(ix0,iy0-1,iz0))
	 rhilpz = difz(iz0)*(sz(ix0,iy0,iz0+1)-sz(ix0,iy0,iz0-1))
	 write(6,48) '   rhilpx: ',rhilpx
	 write(6,48) '   rhilpy: ',rhilpy
	 write(6,48) '   rhilpz: ',rhilpz
      endif
c
      do 90 iz = 2,nzz-1
      do 90 iy = 2,ny1
      ixanf = 2 + mod(iz+iy+igrid,2)
      do 90 ix = ixanf,nx1,2
	  help(ix,iy,iz) = difx(ix)*(sx(ix+1,iy,iz)-sx(ix-1,iy,iz))
     +                     + dify(iy)*(sy(ix,iy+1,iz)-sy(ix,iy-1,iz))
     +                     + difz(iz)*(sz(ix,iy,iz+1)-sz(ix,iy,iz-1))
   90 continue
      if (istep .lt. ivisrho) then
	do 100 iz = 2,nzz-1
        do 100 iy = 2,ny1
	ixanf = 2 + mod(iz+iy+igrid,2)
	do 100 ix = ixanf,nx1,2
	   rho(ix,iy,iz) = rho(ix,iy,iz) - ddt * help(ix,iy,iz)
  100   continue
      else
	do 110 iz = 2,nzz-1
        do 110 iy = 2,ny1
	ixanf = 2 + mod(iz+iy+igrid,2)
	do 110 ix = ixanf,nx1,2
	     rho(ix,iy,iz) = rho(ix,iy,iz) - dt*help(ix,iy,iz)
  110   continue

	do 120 iz = 2,nzz-1
        do 120 iy = 2,ny1
	do 120 ix = 2,nx1
	    feldx(ix,iy,iz) = (rho(ix+1,iy,iz)-rho(ix,iy,iz))
     +                         *(rho(ix,iy,iz)-rho(ix-1,iy,iz))
	    feldy(ix,iy,iz) = (rho(ix,iy+1,iz)-rho(ix,iy,iz))
     +                         *(rho(ix,iy,iz)-rho(ix,iy-1,iz))
	    feldz(ix,iy,iz) = (rho(ix,iy,iz+1)-rho(ix,iy,iz))
     +                         *(rho(ix,iy,iz)-rho(ix,iy,iz-1))
  120   continue

        call fbound0
      
	do 127 iz = 2,nzz-1
	do 127 iy = 2,ny1
	ixanf = 2 + mod(iz+iy+igrid,2)
	do 127 ix = ixanf,nx1,2
	   if( ( feldx(ix,iy,iz).lt.0 .and.
     +         (feldx(ix+1,iy,iz).lt.0 .or. feldx(ix-1,iy,iz).lt.0) )
     +      .or. ( feldy(ix,iy,iz).lt.0 .and.
     +         (feldy(ix,iy+1,iz).lt.0 .or. feldy(ix,iy-1,iz).lt.0) )
     +      .or. ( feldz(ix,iy,iz).lt.0 .and.
     +         (feldz(ix,iy,iz+1).lt.0 .or. feldz(ix,iy,iz-1).lt.0) ) )
     +      help(ix,iy,iz) =  help(ix,iy,iz)
     +         -2./ddt*( -(visx+visy+visz)*rho(ix,iy,iz)
     +         + visx*( meanpx(ix)*rho(ix+1,iy,iz)
     +                                 + meanmx(ix)*rho(ix-1,iy,iz) )
     +         + visy*(meanpy(iy)*rho(ix,iy+1,iz)
     +                                 + meanmy(iy)*rho(ix,iy-1,iz) )
     +         + visz*(meanpz(iz)*rho(ix,iy,iz+1)
     +                                 + meanmz(iz)*rho(ix,iy,iz-1) ) )
  127   continue
	do 130 iz = 2,nzz-1
	do 130 iy = 2,ny1
	ixanf = 2 + mod(iz+iy+igrid,2)
	do 130 ix = ixanf,nx1,2
	    rho(ix,iy,iz) = rho(ix,iy,iz) - dt*help(ix,iy,iz)
  130   continue
      end if
c
      if (lndo) then
	 write(6,48) '     help: ',help(ix0,iy0,iz0)
	 write(6,27)
	 write(6,48) ' rho (after advancement): ',rho(ix0,iy0,iz0)
	 write(6,27)
	 write(6,27)
      endif
c
c************************************************
c    2. momentum equations
c    ------------------------
c************************************************
c
      do 220 iz = 1,nzz
      do 220 iy = 1,ny
      ixanf = 1 + mod(iz+iy+igrid,2)
      do 220 ix = ixanf,nx,2	  
         help(ix,iy,iz) = u(ix,iy,iz)**gamma + 0.5*(
     +                        + bx(ix,iy,iz)*bx(ix,iy,iz)
     +                        + by(ix,iy,iz)*by(ix,iy,iz)
     +                        + bz(ix,iy,iz)*bz(ix,iy,iz) ) 
  220 continue
c
      do 240 iz = 1,nzz
      do 240 iy = 1,ny
      ixanf = 1 + mod(iz+iy+igrid,2)
      do 240 ix = ixanf,nx,2
	 hilfx(ix,iy,iz) = help(ix,iy,iz)
     +                   + sx(ix,iy,iz)*vx(ix,iy,iz)
     +                   - bx(ix,iy,iz)*bx(ix,iy,iz)
	 hilfy(ix,iy,iz) = help(ix,iy,iz)
     +                   + sy(ix,iy,iz)*vy(ix,iy,iz)
     +                   - by(ix,iy,iz)*by(ix,iy,iz)
	 hilfz(ix,iy,iz) = help(ix,iy,iz)
     +                   + sz(ix,iy,iz)*vz(ix,iy,iz)
     +                   - bz(ix,iy,iz)*bz(ix,iy,iz)
	 feldx(ix,iy,iz) = sy(ix,iy,iz)*vz(ix,iy,iz)
     +                    - by(ix,iy,iz)*bz(ix,iy,iz)
	 feldy(ix,iy,iz) = sz(ix,iy,iz)*vx(ix,iy,iz)
     +                    - bz(ix,iy,iz)*bx(ix,iy,iz)
	 feldz(ix,iy,iz) = sx(ix,iy,iz)*vy(ix,iy,iz)
     +                    - bx(ix,iy,iz)*by(ix,iy,iz)
  240 continue
c
c-----------------------------------------------------
c     Initialize artificial frictional force
c-----------------------------------------------------
c 
c      do 57 iz = 1,nz
c         write(*,*) 'muprof:', z(iz), muprof(iz)
c   57 continue
     

      do 245 iz = 2,nzz-1
      do 245 iy = 2,ny1
      ixanf = 2 + mod(iz+iy+igrid,2)
      do 245 ix = ixanf,nx1,2
	 hilfx(ix,iy,iz)=difx(ix)*(hilfx(ix+1,iy,iz)-hilfx(ix-1,iy,iz))
     +                 + dify(iy)*(feldz(ix,iy+1,iz)-feldz(ix,iy-1,iz))
     +                 + difz(iz)*(feldy(ix,iy,iz+1)-feldy(ix,iy,iz-1))
     +                 + muprof(ix,iy,iz)*sx(ix,iy,iz)

	 hilfy(ix,iy,iz)=dify(iy)*(hilfy(ix,iy+1,iz)-hilfy(ix,iy-1,iz))
     +                 + difz(iz)*(feldx(ix,iy,iz+1)-feldx(ix,iy,iz-1))
     +                 + difx(ix)*(feldz(ix+1,iy,iz)-feldz(ix-1,iy,iz))
     +                 + muprof(ix,iy,iz)*sy(ix,iy,iz)

         hilfz(ix,iy,iz)=difz(iz)*(hilfz(ix,iy,iz+1)-hilfz(ix,iy,iz-1))
     +                 + difx(ix)*(feldy(ix+1,iy,iz)-feldy(ix-1,iy,iz))
     +                 + dify(iy)*(feldx(ix,iy+1,iz)-feldx(ix,iy-1,iz))
     +                 + muprof(ix,iy,iz)*sz(ix,iy,iz)

  245 continue
      do 250 iz = 2,nzz-1
      do 250 iy = 2,ny1
      ixanf = 2 + mod(iz+iy+igrid,2)
      do 250 ix = ixanf,nx1,2
	 sx(ix,iy,iz) = sx(ix,iy,iz) - dt*hilfx(ix,iy,iz)
	 sy(ix,iy,iz) = sy(ix,iy,iz) - dt*hilfy(ix,iy,iz)
	 sz(ix,iy,iz) = sz(ix,iy,iz) - dt*hilfz(ix,iy,iz)
  250 continue
c
      do 260 iz = 2,nzz-1
      do 260 iy = 2,ny1
      do 260 ix = 2,nx1
         feldx(ix,iy,iz) = (sx(ix+1,iy,iz)-sx(ix,iy,iz))
     +                     *(sx(ix,iy,iz)-sx(ix-1,iy,iz))
	 feldy(ix,iy,iz) = (sx(ix,iy+1,iz)-sx(ix,iy,iz))
     +                     *(sx(ix,iy,iz)-sx(ix,iy-1,iz))
	 feldz(ix,iy,iz) = (sx(ix,iy,iz+1)-sx(ix,iy,iz))
     +                     *(sx(ix,iy,iz)-sx(ix,iy,iz-1))
  260 continue
  
      call fbound0
      
      do 267 iz = 2,nzz-1
      do 267 iy = 2,ny1
      ixanf = 2 + mod(iz+iy+igrid,2)
      do 267 ix = ixanf,nx1,2
	if( ( feldx(ix,iy,iz).lt.0 .and.
     +          (feldx(ix+1,iy,iz).lt.0 .or. feldx(ix-1,iy,iz).lt.0) )
     +     .or.( feldy(ix,iy,iz).lt.0 .and.
     +          (feldy(ix,iy+1,iz).lt.0 .or. feldy(ix,iy-1,iz).lt.0) )
     +     .or.( feldz(ix,iy,iz).lt.0 .and.
     +          (feldz(ix,iy,iz+1).lt.0 .or. feldz(ix,iy,iz-1).lt.0) ) )
     +   hilfx(ix,iy,iz) = hilfx(ix,iy,iz)
     +        -2./ddt*(  -(visx+visy+visz)*sx(ix,iy,iz)
     +    + visx*(meanpx(ix)*sx(ix+1,iy,iz)+meanmx(ix)*sx(ix-1,iy,iz))
     +    + visy*(meanpy(iy)*sx(ix,iy+1,iz)+meanmy(iy)*sx(ix,iy-1,iz))
     +    + visz*(meanpz(iz)*sx(ix,iy,iz+1)+meanmz(iz)*sx(ix,iy,iz-1)) )
  267 continue
c
      do 360 iz = 2,nzz-1
      do 360 iy = 2,ny1
      do 360 ix = 2,nx1
	 feldx(ix,iy,iz) = (sy(ix+1,iy,iz)-sy(ix,iy,iz))
     +                     *(sy(ix,iy,iz)-sy(ix-1,iy,iz))
	 feldy(ix,iy,iz) = (sy(ix,iy+1,iz)-sy(ix,iy,iz))
     +                     *(sy(ix,iy,iz)-sy(ix,iy-1,iz))
	 feldz(ix,iy,iz) = (sy(ix,iy,iz+1)-sy(ix,iy,iz))
     +                     *(sy(ix,iy,iz)-sy(ix,iy,iz-1))
  360 continue
  
      call fbound0
      
      do 367 iz = 2,nzz-1
      do 367 iy = 2,ny1
      ixanf = 2 + mod(iz+iy+igrid,2)
      do 367 ix = ixanf,nx1,2
	 if( ( feldx(ix,iy,iz).lt.0 .and.
     +          (feldx(ix+1,iy,iz).lt.0 .or. feldx(ix-1,iy,iz).lt.0) )
     +     .or.( feldy(ix,iy,iz).lt.0 .and.
     +          (feldy(ix,iy+1,iz).lt.0 .or. feldy(ix,iy-1,iz).lt.0) )
     +     .or.( feldz(ix,iy,iz).lt.0 .and.
     +          (feldz(ix,iy,iz+1).lt.0 .or. feldz(ix,iy,iz-1).lt.0) ) )
     +   hilfy(ix,iy,iz) =  hilfy(ix,iy,iz)
     +        -2./ddt*(  -(visx+visy+visz)*sy(ix,iy,iz)
     +    + visx*(meanpx(ix)*sy(ix+1,iy,iz)+meanmx(ix)*sy(ix-1,iy,iz))
     +    + visy*(meanpy(iy)*sy(ix,iy+1,iz)+meanmy(iy)*sy(ix,iy-1,iz))
     +    + visz*(meanpz(iz)*sy(ix,iy,iz+1)+meanmz(iz)*sy(ix,iy,iz-1)) )
  367 continue
      do 460 iz = 2,nzz-1
      do 460 iy = 2,ny1
      do 460 ix = 2,nx1
	 feldx(ix,iy,iz) = (sz(ix+1,iy,iz)-sz(ix,iy,iz))
     +                     *(sz(ix,iy,iz)-sz(ix-1,iy,iz))
	 feldy(ix,iy,iz) = (sz(ix,iy+1,iz)-sz(ix,iy,iz))
     +                     *(sz(ix,iy,iz)-sz(ix,iy-1,iz))
	 feldz(ix,iy,iz) = (sz(ix,iy,iz+1)-sz(ix,iy,iz))
     +                     *(sz(ix,iy,iz)-sz(ix,iy,iz-1))
  460 continue
  
      call fbound0
      
      do 467 iz = 2,nzz-1
      do 467 iy = 2,ny1
      ixanf = 2 + mod(iz+iy+igrid,2)
      do 467 ix = ixanf,nx1,2
	if( ( feldx(ix,iy,iz).lt.0 .and.
     +          (feldx(ix+1,iy,iz).lt.0 .or. feldx(ix-1,iy,iz).lt.0) )
     +     .or.( feldy(ix,iy,iz).lt.0 .and.
     +          (feldy(ix,iy+1,iz).lt.0 .or. feldy(ix,iy-1,iz).lt.0) )
     +     .or.( feldz(ix,iy,iz).lt.0 .and.
     +          (feldz(ix,iy,iz+1).lt.0 .or. feldz(ix,iy,iz-1).lt.0) ) )
     +   hilfz(ix,iy,iz) =  hilfz(ix,iy,iz)
     +        -2./ddt*(  -(visx+visy+visz)*sz(ix,iy,iz)
     +    + visx*(meanpx(ix)*sz(ix+1,iy,iz)+meanmx(ix)*sz(ix-1,iy,iz))
     +    + visy*(meanpy(iy)*sz(ix,iy+1,iz)+meanmy(iy)*sz(ix,iy-1,iz))
     +    + visz*(meanpz(iz)*sz(ix,iy,iz+1)+meanmz(iz)*sz(ix,iy,iz-1)) )
  467 continue
      do 500 iz = 2,nzz-1
      do 500 iy = 2,ny1
      ixanf = 2 + mod(iz+iy+igrid,2)
      do 500 ix = ixanf,nx1,2
	 sx(ix,iy,iz) = sx(ix,iy,iz) - dt*hilfx(ix,iy,iz)
	 sy(ix,iy,iz) = sy(ix,iy,iz) - dt*hilfy(ix,iy,iz)
	 sz(ix,iy,iz) = sz(ix,iy,iz) - dt*hilfz(ix,iy,iz)
  500 continue
c
c************************************************
c    3.a resistivity
c    -----------------
c************************************************
c
        do 499 iz = 1,nzz
	do 499 iy = 1,ny
	do 499 ix = 1,nx
	   feldx(ix,iy,iz) = 0.
	   feldy(ix,iy,iz) = 0.
	   feldz(ix,iy,iz) = 0.
	   hilf(ix,iy,iz) = 0.
  499   continue

        do 501 iz = 2,nzz-1
        do 501 iy = 2,ny1
	ixanf = 2 + mod(iz+iy+igrid,2)
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
c
c************************************************
c    3. magnetic field
c    -------------------
c************************************************
c
c     hilf = v x b
c
      do 510 iz = 1,nzz
      do 510 iy = 1,ny 
      do 510 ix = 1,nx
	 hilfx(ix,iy,iz) = 0.
	 hilfy(ix,iy,iz) = 0.
	 hilfz(ix,iy,iz) = 0.
  510 continue
c
      do 520 iz = 1,nzz
      do 520 iy = 1,ny
      ixanf = 1 + mod(iz+iy+igrid,2)
      do 520 ix = ixanf,nx,2
	 hilfx(ix,iy,iz) = vy(ix,iy,iz)*bz(ix,iy,iz)
     +                     - vz(ix,iy,iz)*by(ix,iy,iz) 
	 hilfy(ix,iy,iz) = vz(ix,iy,iz)*bx(ix,iy,iz)
     +                     - vx(ix,iy,iz)*bz(ix,iy,iz) 
	 hilfz(ix,iy,iz) = vx(ix,iy,iz)*by(ix,iy,iz)
     +                     - vy(ix,iy,iz)*bx(ix,iy,iz) 
  520 continue
c
c      feld = rot b
c      help notwendig fuer implizite berechnung von b
c
      do 540 iz = 2,nzz-1
      do 540 iy = 2,ny1 
      ixanf = 2 + mod(iz+iy+igrid,2)
      do 540 ix = ixanf,nx1,2
	 help(ix,iy,iz) = res(ix,iy,iz)*(
     +                        ddifx(ix) + ddify(iy) + ddifz(iz))
  540 continue
c
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
         write(6,*) 'Subroutine leap'
	 write(6,*) 'Magnetic field diagnostics'
         write(6,47) '    istep: ',istep
         write(6,48) '     time: ',time
         write(6,48) 'time step: ',dt
         write(6,48) 'double dt: ',ddt
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
c        term3: eta Laplacian b
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
         term3(1) =
     +        + res(ix0,iy0,iz0)*(
     +              (ddifpx(ix0)*bx(ix0+1,iy0,iz0)+
     +               ddifmx(ix0)*bx(ix0-1,iy0,iz0))   
     +            + (ddifpy(iy0)*bx(ix0,iy0+1,iz0)+
     +               ddifmy(iy0)*bx(ix0,iy0-1,iz0))  
     +            + (ddifpz(iz0)*bx(ix0,iy0,iz0+1)+
     +               ddifmz(iz0)*bx(ix0,iy0,iz0-1)) )
     +        - help(ix0,iy0,iz0)*bx(ix0,iy0,iz0)
c
         term1(2) = 
     +             difz(iz0)*( 
     +                  hilfx(ix0,iy0,iz0+1)-hilfx(ix0,iy0,iz0-1) ) 
     +           - difx(ix0)*( 
     +                  hilfz(ix0+1,iy0,iz0)-hilfz(ix0-1,iy0,iz0) )
         term2(2) =
     +        - feldx(ix0,iy0,iz0)*difz(iz0)*(
     +                        res(ix0,iy0,iz0+1)-res(ix0,iy0,iz0-1))
     +        + feldz(ix0,iy0,iz0)*difx(ix0)*(
     +                        res(ix0+1,iy0,iz0)-res(ix0-1,iy0,iz0))
         term3(2) = 
     +        + res(ix0,iy0,iz0)*(
     +              (ddifpx(ix0)*by(ix0+1,iy0,iz0)+
     +               ddifmx(ix0)*by(ix0-1,iy0,iz0))
     +            + (ddifpy(iy0)*by(ix0,iy0+1,iz0)+
     +               ddifmy(iy0)*by(ix0,iy0-1,iz0))
     +            + (ddifpz(iz0)*by(ix0,iy0,iz0+1)+
     +               ddifmz(iz0)*by(ix0,iy0,iz0-1)) )
     +        - help(ix0,iy0,iz0)*by(ix0,iy0,iz0)
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
         term3(3) =
     +        + res(ix0,iy0,iz0)*(
     +              (ddifpx(ix0)*bz(ix0+1,iy0,iz0)+
     +               ddifmx(ix0)*bz(ix0-1,iy0,iz0))
     +            + (ddifpy(iy0)*bz(ix0,iy0+1,iz0)+
     +               ddifmy(iy0)*bz(ix0,iy0-1,iz0))
     +            + (ddifpz(iz0)*bz(ix0,iy0,iz0+1)+
     +               ddifmz(iz0)*bz(ix0,iy0,iz0-1)) )
     +        - help(ix0,iy0,iz0)*bz(ix0,iy0,iz0)
c
         dterm = 1 + ddt*help(ix0,iy0,iz0)
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
         write(6,48) 'eta Laplacian b:',term3(1)
         write(6,27)
         write(6,27)
         write(6,*) 'by'
         write(6,*) '--'
         write(6,48) 'curl ( v X b ):',term1(2)
         write(6,48) '- grad ( eta ) X j:',term2(2)
         write(6,48) 'eta Laplacian b:',term3(2)
         write(6,27)
         write(6,27)
         write(6,*) 'bz'
         write(6,*) '--'
         write(6,48) 'curl ( v X b ):',term1(3)
         write(6,48) '- grad ( eta ) X j:',term2(3)
         write(6,48) 'eta Laplacian b:',term3(3)
         write(6,27)
         write(6,27)
         write(6,27)
         write(6,48) 'Divisor term:',dterm
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
      do 560 iz = 2,nzz-1
      do 560 iy = 2,ny1
      ixanf = 2 + mod(iz+iy+igrid,2)
      do 560 ix = ixanf,nx1,2
	bx(ix,iy,iz) = ( bx(ix,iy,iz) + ddt*(
     +           dify(iy)*( hilfz(ix,iy+1,iz)-hilfz(ix,iy-1,iz) )
     +         - difz(iz)*( hilfy(ix,iy,iz+1)-hilfy(ix,iy,iz-1) )
     +     - feldz(ix,iy,iz)*dify(iy)*(res(ix,iy+1,iz)-res(ix,iy-1,iz))     
     +     + feldy(ix,iy,iz)*difz(iz)*(res(ix,iy,iz+1)-res(ix,iy,iz-1))
     +     + res(ix,iy,iz)*(
     +           (ddifpx(ix)*bx(ix+1,iy,iz)+ddifmx(ix)*bx(ix-1,iy,iz))   
     +         + (ddifpy(iy)*bx(ix,iy+1,iz)+ddifmy(iy)*bx(ix,iy-1,iz))  
     +         + (ddifpz(iz)*bx(ix,iy,iz+1)+ddifmz(iz)*bx(ix,iy,iz-1)) )
     +     - help(ix,iy,iz)*bx(ix,iy,iz)
     +                                  ) ) /(1 + ddt*help(ix,iy,iz))
	by(ix,iy,iz) = ( by(ix,iy,iz) + ddt*(
     +          difz(iz)*( hilfx(ix,iy,iz+1)-hilfx(ix,iy,iz-1) ) 
     +        - difx(ix)*( hilfz(ix+1,iy,iz)-hilfz(ix-1,iy,iz) )
     +     - feldx(ix,iy,iz)*difz(iz)*(res(ix,iy,iz+1)-res(ix,iy,iz-1))
     +     + feldz(ix,iy,iz)*difx(ix)*(res(ix+1,iy,iz)-res(ix-1,iy,iz))
     +     + res(ix,iy,iz)*(
     +           (ddifpx(ix)*by(ix+1,iy,iz)+ddifmx(ix)*by(ix-1,iy,iz))
     +         + (ddifpy(iy)*by(ix,iy+1,iz)+ddifmy(iy)*by(ix,iy-1,iz))
     +         + (ddifpz(iz)*by(ix,iy,iz+1)+ddifmz(iz)*by(ix,iy,iz-1)) )
     +     - help(ix,iy,iz)*by(ix,iy,iz)
     +                                  ) ) / (1 + ddt*help(ix,iy,iz))
	bz(ix,iy,iz) = ( bz(ix,iy,iz) + ddt*(
     +          difx(ix)*( hilfy(ix+1,iy,iz)-hilfy(ix-1,iy,iz) )
     +        - dify(iy)*( hilfx(ix,iy+1,iz)-hilfx(ix,iy-1,iz) )
     +     - feldy(ix,iy,iz)*difx(ix)*(res(ix+1,iy,iz)-res(ix-1,iy,iz))
     +     + feldx(ix,iy,iz)*dify(iy)*(res(ix,iy+1,iz)-res(ix,iy-1,iz))
     +     + res(ix,iy,iz)*(
     +           (ddifpx(ix)*bz(ix+1,iy,iz)+ddifmx(ix)*bz(ix-1,iy,iz))
     +         + (ddifpy(iy)*bz(ix,iy+1,iz)+ddifmy(iy)*bz(ix,iy-1,iz))
     +         + (ddifpz(iz)*bz(ix,iy,iz+1)+ddifmz(iz)*bz(ix,iy,iz-1)) )
     +     - help(ix,iy,iz)*bz(ix,iy,iz)
     +                                  ) ) / (1 + ddt*help(ix,iy,iz))
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
      do 820 iz = 1,nzz
      do 820 iy = 1,ny
      ixanf = 1 + mod(iz+iy+igrid,2)
      do 820 ix = ixanf,nx,2
	 hilfx(ix,iy,iz) = u(ix,iy,iz)*vx(ix,iy,iz)
	 hilfy(ix,iy,iz) = u(ix,iy,iz)*vy(ix,iy,iz)
	 hilfz(ix,iy,iz) = u(ix,iy,iz)*vz(ix,iy,iz)
  820 continue
      chelp = eingam/gamma/3.0**eingam
      do 824 iz = 2,nzz-1
      do 824 iy = 2,ny1
      do 824 ix = 2,nx
	 hilpx(ix,iy,iz) = 0.0
	 hilpy(ix,iy,iz) = 0.0
	 hilpz(ix,iy,iz) = 0.0
 824  continue
      do 825 iz = 2,nzz-1
      do 825 iy = 2,ny1
      ixanf = 2 + mod(iz+iy+igrid,2)
      do 825 ix = ixanf,nx1,2
	 hilpx(ix,iy,iz) = difx(ix)*(hilfx(ix+1,iy,iz)-hilfx(ix-1,iy,iz))
	 hilpy(ix,iy,iz) = dify(iy)*(hilfy(ix,iy+1,iz)-hilfy(ix,iy-1,iz))
	 hilpz(ix,iy,iz) = difz(iz)*(hilfz(ix,iy,iz+1)-hilfz(ix,iy,iz-1))
	 help(ix,iy,iz) = chelp*res(ix,iy,iz)*hilf(ix,iy,iz)*(
     +             ( meanpx(ix)*u(ix+1,iy,iz)+meanmx(ix)*u(ix-1,iy,iz) )
     +            +( meanpy(iy)*u(ix,iy+1,iz)+meanmy(iy)*u(ix,iy-1,iz) )
     +            +( meanpz(iz)*u(ix,iy,iz+1)+meanmz(iz)*u(ix,iy,iz-1) )
     +                           )**eingam 
 825  continue
      do 826 iz = 2,nzz-1
      do 826 iy = 2,ny1
      ixanf = 2 + mod(iz+iy+igrid,2)
      do 826 ix = ixanf,nx1,2
	 help(ix,iy,iz) = help(ix,iy,iz) + 
     +                 hilpx(ix,iy,iz)+hilpy(ix,iy,iz)+hilpz(ix,iy,iz)
 826  continue
c
c
c     Additional pressure diagnostics
      if (lndo) then
	 write(6,27)
	 write(6,27)
	 write(6,*) 'Subroutine leap'
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
	 write(6,*) 'Derived quantities'
	 write(6,48) ' help: ',hlpstore
	 write(6,48) 'hilpx: ',hilpx(ix0,iy0,iz0)
	 write(6,48) 'hilpy: ',hilpy(ix0,iy0,iz0)
	 write(6,48) 'hilpz: ',hilpz(ix0,iy0,iz0)
	 write(6,48) ' hilf: ',hilf(ix0,iy0,iz0)
	 write(6,27)
	 prezok = hilpx(ix0,iy0,iz0) + hilpy(ix0,iy0,iz0) + 
     +            hilpz(ix0,iy0,iz0)
	 write(6,48) ' Sum of hilp terms: ',prezok
	 write(6,48) '  Advancement term: ',ddt*prezok
	 write(6,48) 'Estimate for new u: ',u(ix0,iy0,iz0) - ddt*prezok
	 write(6,27)
      endif

c

      if (istep .lt. ivisu) then
	do 900 iz = 2,nzz-1
	do 900 iy = 2,ny1
	ixanf = 2 + mod(iz+iy+igrid,2)
	do 900 ix = ixanf,nx1,2
	   u(ix,iy,iz) = u(ix,iy,iz) - ddt * help(ix,iy,iz)
  900   continue
      else
	do 910 iz = 2,nzz-1
	do 910 iy = 2,ny1
	ixanf = 2 + mod(iz+iy+igrid,2)
	do 910 ix = ixanf,nx1,2
	   u(ix,iy,iz) = u(ix,iy,iz) - dt*help(ix,iy,iz)
  910   continue
c
	do 920 iz = 2,nzz-1
	do 920 iy = 2,ny1
	do 920 ix = 2,nx1
	   feldx(ix,iy,iz) = (u(ix+1,iy,iz)-u(ix,iy,iz))
     +                         *(u(ix,iy,iz)-u(ix-1,iy,iz))
	   feldy(ix,iy,iz) = (u(ix,iy+1,iz)-u(ix,iy,iz))
     +                         *(u(ix,iy,iz)-u(ix,iy-1,iz))
	   feldz(ix,iy,iz) = (u(ix,iy,iz+1)-u(ix,iy,iz))
     +                         *(u(ix,iy,iz)-u(ix,iy,iz-1))
  920   continue
  
         call fbound0
      
	 do 927 iz = 2,nzz-1
	 do 927 iy = 2,ny1
	 ixanf = 2 + mod(iz+iy+igrid,2)
	 do 927 ix = ixanf,nx1,2
	  if( ( feldx(ix,iy,iz).lt.0 .and.
     +          (feldx(ix+1,iy,iz).lt.0 .or. feldx(ix-1,iy,iz).lt.0) )
     +       .or.( feldy(ix,iy,iz).lt.0 .and.
     +          (feldy(ix,iy+1,iz).lt.0 .or. feldy(ix,iy-1,iz).lt.0) )
     +       .or.( feldz(ix,iy,iz).lt.0 .and.
     +          (feldz(ix,iy,iz+1).lt.0 .or. feldz(ix,iy,iz-1).lt.0) ) )
     +     help(ix,iy,iz) =  help(ix,iy,iz)
     +         -2./ddt*( -(visx+visy+visz)*u(ix,iy,iz)
     +      + visx*(meanpx(ix)*u(ix+1,iy,iz)+meanmx(ix)*u(ix-1,iy,iz))
     +      + visy*(meanpy(iy)*u(ix,iy+1,iz)+meanmy(iy)*u(ix,iy-1,iz))
     +      + visz*(meanpz(iz)*u(ix,iy,iz+1)+meanmz(iz)*u(ix,iy,iz-1)) )
  927    continue
	 do 930 iz = 2,nzz-1
	 do 930 iy = 2,ny1
	   ixanf = 2 + mod(iz+iy+igrid,2)
	   do 930 ix = ixanf,nx1,2
	     u(ix,iy,iz) = u(ix,iy,iz) - dt*help(ix,iy,iz)
  930    continue
      end if
      if (mod((istep+1), 100) .eq. 0) then
        write(*,*) 'Pressure, time:', time
        do 935 ix = nx,nx-11,-1
  	  write(*,45) x(ix), u(ix,2,2), 
     +              hilpx(ix,2,2),hilpy(ix,2,2),hilpz(ix,2,2)
 935    continue
      endif
c
c-----------------------------------------------------------------
c     Activate relaxation switch if necessary.
      pinf=1.5
      if ((istep+1) .eq. rsw) then
        nrelax=newrelax
        pinf=4.
      endif
      if ((istep+1) .eq. rsw1) then
        nrelax=newrelax1
        pinf=4.
      endif
      if ((istep+1) .eq. rsw2) then
        nrelax=newrelax2
        pinf=4.
      endif
      if ((istep+1) .eq. rsw3) nrelax=newrelax3
      if ((istep+1) .eq. rsw4) nrelax=newrelax4
      if ((istep+1) .eq. rsw5) nrelax=newrelax5
      if ((istep+1) .eq. rsw6) nrelax=newrelax6
      if ((istep+1) .eq. rsw7) nrelax=newrelax7
      if ((istep+1) .eq. rsw8) nrelax=newrelax8
c
c     Apply `ballistic relaxation.'
c
      if ((mod((istep+1), nrelax) .eq. 0) .and. relax  
     +        .and.   (istep+1 .le. nmrelax) )  then
	 do 937 iz = 1,nz
	 do 937 iy = 1,ny
	 do 937 ix = 1,nx
	   sx(ix,iy,iz) = 0.0
	   sy(ix,iy,iz) = 0.0
	   sz(ix,iy,iz) = 0.0
 937	continue
        do 940 ix = 1,nx
        do 940 iy = 1,ny
        do 940 iz = 1,nz
	  rho(ix,iy,iz)=2.*u(ix,iy,iz)**gamma
 940    continue
        pmin=100.
        do 948 ix = 1,nx
        do 948 iy = 1,ny
        do 948 iz = 1,nz
          if (rho(ix,iy,iz) .lt. pmin) then
            pmin=rho(ix,iy,iz)
            ixx=ix
            iyy=iy
            izz=iz
          endif
 948    continue
c
        do 949 ix = 1,nx
        do 949 iy = 1,ny
        do 949 iz = 1,nz
          rho(ix,iy,iz)=rho(ix,iy,iz)-pmin+pinf
 949    continue
        do 950 ix = 1,nx
        do 950 iy = 1,ny
        do 950 iz = 1,nz
          u(ix,iy,iz) = (rho(ix,iy,iz)/2.0)**gaminv
 950    continue


	write(6,27)
	write(6,*) 'Ballistic relaxation applied'
	write(6,47) 'Current value of istep = ',istep
	write(6,48) 'Current value of  time = ',time
	write(6,47) 'Relaxation effects evident at istep = ',
     +               istep+1
	write(6,48) 'Projected time = ',(istep+1)*dt
	write(6,27)
      endif
c
c------------------------------------------------------------------
c

  27  format((/))
  45  format(1x,f7.2, 4f10.4)
  47  format( 1x,A,i5)
  48  format( 1x,A,f14.5)
  50  format(1x,'Indices: (',i3,',',i3,',',i3,')')
  55  format(1x,'Location: (',f8.3,',',f8.3,',',f8.3,')')


      return
      end
