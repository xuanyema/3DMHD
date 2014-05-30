	subroutine initcon
c***********************************************************************
	include 'misflin'
c
	integer  ix,iy,iz
	real     gaminv,phirad,psirad,eta0,sinpsi,cospsi,
     +           pprof(nx),n1,n2,ylim,zlim,ytrans,ztrans,
     +           mu0,zmu,dzmu
c ...................................................
c  initial configuration
c ...................................................
      eta0   = 0.0
      if (etasw.eq.0)  eta0 = eta
      gaminv = 1./gamma
      phirad = 2.*pi*phi/360.
      psirad = 2.*pi*psi/360.
      delbz  = 0.5*sqrt( 1.+bmsh**2-2.*bmsh*cos(phirad) )
      by0    = 0.5*bmsh*sin(phirad)/delbz
      bz0    = 0.25*(1. - bmsh**2)/delbz
      p0     = pmsp + 0.5*(1.0-bmsh**2)
      bzmsp  = bz0 + delbz
      bzmsh  = bz0 - delbz
      pmsh   = pmsp + 1.0 -bmsh**2
      betasp = pmsp
      betash = (pmsp+1.0-bmsh**2)/bmsh**2
      vamsp  = sqrt( 1.0/(rho0-delrho) )
      vamsh  = sqrt( bmsh**2/(rho0+delrho) )
      csmsp  = sqrt( gamma/2.0*pmsp/(rho0-delrho) )
      csmsh  = sqrt( gamma/2.0*(pmsp+1.0-bmsh**2)/(rho0+delrho) )
c
c  equilibrium profiles
c   help used auxiliary to store pressure
      do 50 ix = 1,nx
	 rhoprof(ix) = rho0 + delrho*tanh((x(ix)-xrho)/dxrho)
         pprof(ix) = p0 + 2.*bz0*delbz*tanh(x(ix))
     +                      + (1.-kappa)*delbz*delbz
     +                          /cosh(x(ix))/cosh(x(ix))
	 bxprof(ix)  = 0.0
	 byprof(ix)  = sqrt(by0*by0+kappa*delbz*delbz
     +                                /cosh(x(ix))/cosh(x(ix)) )
	 bzprof(ix)  = bz0 - delbz*tanh(x(ix))
	 sxprof(ix)  = 0.0
	 syprof(ix)  = vy0*rhoprof(ix)*tanh(x(ix))
	 szprof(ix)  = 0.0
   50 continue
      write(*,*) ' In Anfang:', rhoprof(2), bzprof(2), byprof(2)
      do 55 ix = 1,nx
	 uprof(ix) = (pprof(ix)/2.)**gaminv
   55 continue
      mu0 = 1.
      dzmu = 3.
      zmu = 3.*zmax/4.
      do 57 iz = 1,nz
	 muprof(iz)  = mu0*
     +       *0.5*(-tanh((z(iz)+zmu)/dzmu)+tanh((z(iz)-zmu)/dzmu)+2.)
         write(*,*) 'muprof:', z(iz), muprof(iz)
   57 continue
      
      if (start .eqv. .false.) return 

      do 60 iz = 1,nz
	do 60 iy = 1,ny
	do 60 ix = 1,nx
	  rho(ix,iy,iz) = rhoprof(ix)
	  help(ix,iy,iz)= pprof(ix)
	  bx(ix,iy,iz)  = bxprof(ix)
c     +               +vx1*pi/ymax*cos(pi*x(ix)/2./xmax)
c     +                              *cos(pi*z(iz)/zmax)
	  by(ix,iy,iz)  = byprof(ix)
	  bz(ix,iy,iz)  = bzprof(ix)
c     +           +vx1*pi/2./xmax*sin(pi*x(ix)/2./xmax)
c     +                 *sin(pi*z(iz)/zmax)
	  sx(ix,iy,iz)  = sxprof(ix)
	  sy(ix,iy,iz)  = syprof(ix)
	  sz(ix,iy,iz)  = szprof(ix)
	  prof(ix,iy,iz) = 1.0
c	  prof(ix,iy,iz) = (1. - tanh(x(ix))**2)
c     +                    *(1. - tanh((z(iz)-0.)/3.)**2)
cc     +       *0.5*( tanh((y(iy)+10.)/4.) - tanh((y(iy)-10.)/4.) )
c     +                    *(1. - tanh(y(iy)/5.)**2)
   60 continue
	write(*,*) ' In Anfanga:', rhoprof(2), sx(2,2,2), bx(2,2,2),vx1
	write(*,*) ' In Anfangb:', pi,xmax,ymax,zmax,x(2),z(2)
        do 499 iz = 1,nz
	do 499 iy = 1,ny
	do 499 ix = 1,nx
	   feldx(ix,iy,iz) = 0.
	   feldy(ix,iy,iz) = 0.
	   feldz(ix,iy,iz) = 0.
  499   continue
	write(*,*) ' In Anfang1:', rhoprof(2), sx(2,2,2)

        n1 = 1.0
        n2 = 4.0
        ylim = 10.
        ytrans = 3.
        zlim =  15.
        ztrans = 3.
        do 502 iz = 2,nz-1
        do 502 iy = 2,ny-1
	do 502 ix = 2,nx-1
	  sx(ix,iy,iz) = sx(ix,iy,iz)
     +             + vx1*n1*pi/ymax * sin(n1*pi/ymax*y(iy))
     +                          /cosh(x(ix)/2.)**2
c     +       *0.5*(tanh((y(iy)+ylim)/ytrans)-tanh((y(iy)-ylim))/ytrans)
     +       *0.5*(tanh((z(iz)+zlim)/ztrans)-tanh((z(iz)-zlim))/ztrans)
     +                          * rho(ix,iy,iz)
	  sy(ix,iy,iz) = sy(ix,iy,iz)
     +             - vx1* cos(n1*pi/ymax*y(iy))
     +                 * tanh(x(ix)/2.)/cosh(x(ix)/2.)**2
c     +       *0.5*(tanh((y(iy)+ylim)/ytrans)-tanh((y(iy)-ylim))/ytrans)
     +       *0.5*(tanh((z(iz)+zlim)/ztrans)-tanh((z(iz)-zlim))/ztrans)
     +                 * rho(ix,iy,iz)
	  sz(ix,iy,iz) = sz(ix,iy,iz)
  502   continue

c
c  perturbations + pressure + resistivity
c   
      do 100 iz = 1,nz
	do 100 iy = 1,ny
	do 100 ix = 1,nx
c  here help is auxiliary for pressure
c	    help(ix,iy,iz) = help(ix,iy,iz) 
c     +      + 0.75/cosh((x(ix)-5.)/2.)**2
c     +           *0.5*( 1.-tanh((z(iz)-0.4*zmax)/3.) )
c     +           *0.5*( tanh((y(iz)+0.75*ymax)/3.) 
c     +                    -tanh((y(iz)-0.75*ymax)/3.) )
c     +           *(1.1/cosh((z(iz)-8.)/2)/cosh((y(iy)-10.)/3)   
c     +              + 1.1/cosh((z(iz)+8.)/2)/cosh((y(iy)+10.)/2)
c     +              + 1./cosh((z(iz)-18.)/2)/cosh((y(iy)+20.)/2)
c     +              + 0.9/cosh((z(iz)-25.)/2)/cosh((y(iy)-5.)/2)
c     +              + 1.2/cosh((z(iz)-40.)/2)/cosh((y(iy)-20.)/2)
c     +              + 0.9/cosh((z(iz)-45.)/2)/cosh((y(iy)+15.)/2) )
	    u(ix,iy,iz)   = (help(ix,iy,iz)/2.)**gaminv
c	    rho(ix,iy,iz) = rho(ix,iy,iz)
c     +      + 0.375/cosh((x(ix)-5.)/2.)**2
c     +           *0.5*( 1.-tanh((z(iz)-0.5*zmax)/3.) )
c     +           *0.5*( tanh((y(iz)+0.75*ymax)/3.)
c     +                    -tanh((y(iz)-0.75*ymax)/3.) )
c     +           *(1.1/cosh((z(iz)-8.)/2)/cosh((y(iy)-10.)/2)
c     +              + 1.1/cosh((z(iz)+8.)/2)/cosh((y(iy)+10.)/2)
c     +              + 1./cosh((z(iz)-18.)/2)/cosh((y(iy)+20.)/2)
c     +              + 0.9/cosh((z(iz)-25.)/2)/cosh((y(iy)-5.)/2)
c     +              + 1.2/cosh((z(iz)-40.)/2)/cosh((y(iy)-20.)/2)
c     +              + 0.9/cosh((z(iz)-45.)/2)/cosh((y(iy)+15.)/2) )
	    res(ix,iy,iz) = eta0*prof(ix,iy,iz)+0.0005
c----------- Init Con's	for bx, sx  (parameters in magin)
c Init Con 1: 
c      etasw=1
c      eta=0.1
c      iend=2000
c      isafe=2000
c      zentr(3)=.true.        for z
c      perio(3)=.false.       for z 
c      lsym(*,3)=.true.       for z
c
c	    bx(ix,iy,iz)  = 0.0
c	    sx(ix,iy,iz)  = 0.0
c Init Con 2: 
c      etasw=2
c      eta=0.04
c      iend=2000
c      isafe=2000
c      zentr(3)=.true.        for z
c      perio(3)=.false.       for z 
c      lsym(*,3)=.true.       for z
c
c	    bx(ix,iy,iz)  = - 0.05*tanh(z(iz))/cosh(0.2*z(iz))
c     +                        /cosh(0.2*y(iy))	    
c	    sx(ix,iy,iz)  = 0.0
c Init Con 3:
c      etasw=2
c      eta=0.04
c      iend=2000
c      isafe=2000
c      zentr(3)=.true.        for z
c      perio(3)=.false.       for z (lsym=.true. for z)
c      lsym(*,3)=.true.       for z
c 
c	    bx(ix,iy,iz)  = 0.0
c	    sx(ix,iy,iz)  = 0.8/cosh(0.3*x(ix))**2
c     +                   /cosh(0.15*z(iz))**2/cosh(0.15*y(iy))**2    
c Init Con 4: 
c      etasw=2
c      eta=0.1
c      iend=2400
c      isafe=2400
c      zentr(3)=.false.        for z
c      perio(3)=.true.         for z (lsym=.false. for z)
c      lsym(*,3)=.false.       for z
c
c	    bx(ix,iy,iz)  = 0.0
c	    sx(ix,iy,iz)  = vx1*pi/ymax*cos(pi/ymax*y(iy))
c     +                      /cosh(x(ix))**2*rho(ix,iy,iz)
c	    sy(ix,iy,iz)  = vy0*tanh(x(ix))*rho(ix,iy,iz)
c     +                   - vx1*sin(pi/ymax*y(iy))
c     +         *2.0*tanh(x(ix))/cosh(x(ix))**2*rho(ix,iy,iz)
  100 continue
c
c    rotate y, z components by psi
c      (recommended for certain wave vectors of perturbation etc.)
      if (psi .ne. 0.0) then 
       sinpsi=sin(psirad)
       cospsi=cos(psirad)
       do 110 iz = 1,nz
       do 110 iy = 1,ny
	do 110 ix = 1,nx
c	    hilfy(ix,iy,iz) = sy(ix,iy,iz)
c	    hilfz(ix,iy,iz) = sz(ix,iy,iz)
	    feldy(ix,iy,iz) = by(ix,iy,iz)
	    feldz(ix,iy,iz) = bz(ix,iy,iz)
  110  continue
       do 120 iz = 1,nz
       do 120 iy = 1,ny
        do 120 ix = 1,nx
c	 sy(ix,iy,iz) = cospsi*hilfy(ix,iy,iz) + sinpsi*hilfz(ix,iy,iz)
c	 sz(ix,iy,iz) = -sinpsi*hilfy(ix,iy,iz) + cospsi*hilfz(ix,iy,iz)
	 by(ix,iy,iz) = cospsi*feldy(ix,iy,iz) + sinpsi*feldz(ix,iy,iz)
	 bz(ix,iy,iz) = -sinpsi*feldy(ix,iy,iz) + cospsi*feldz(ix,iy,iz)
  120  continue
      endif
c
c    boundary values at xmin and xmax
c
	write(*,*) ' In Anfang:', rhoprof(2), sx(2,2,2)
      rhobd(1) = rhoprof(2)
      ubd(1)   = uprof(2)
      pbd(1)   = 2.*uprof(2)**gamma
      vxbd(1)  = sx(2,2,2)/rhoprof(2)
      vybd(1)  = sy(2,2,2)/rhoprof(2)
      vzbd(1)  = sz(2,2,2)/rhoprof(2)
      bxbd(1)  = bx(2,2,2)
      bybd(1)  = by(2,2,2)
      bzbd(1)  = bz(2,2,2)
      rhobd(2) = rhoprof(nx1)
      ubd(2)   = uprof(nx1)
      pbd(2)   = 2.*uprof(nx1)**gamma
      vxbd(2)  = sx(nx1,2,2)/rhoprof(nx1)
      vybd(2)  = sy(nx1,2,2)/rhoprof(nx1)
      vzbd(2)  = sz(nx1,2,2)/rhoprof(nx1)
      bxbd(2)  = bx(nx1,2,2)
      bybd(2)  = by(nx1,2,2)
      bzbd(2)  = bz(nx1,2,2)
      
      write(*,212)  rhobd(1),pbd(1),vxbd(1),vybd(1),vzbd(1),
     +                bxbd(1),bybd(1),bzbd(1),
     +                rhobd(2),pbd(2),vxbd(2),vybd(2),vzbd(2),
     +                bxbd(2),bybd(2),bzbd(2)
  212 format(/1x,'Initial asymptotic values at xmin (1. row) ',
     +           'and xmax (2. row) for: ',
     +       /1x,'    rho      p        vx       vy       vz    ',
     +           '   bx       by       bz    ',
     +       /,8f9.3/,8f9.3)
      write(*,*) 'psi:',psi
      write(*,*) 'MSP: v*k,va*k,cs:',vybd(1),bybd(1)/sqrt(rhobd(1)),
     +            sqrt(0.5*gamma*pbd(1)/rhobd(1))
      write(*,*) 'MSH: v*k,va*k,cs:',vybd(2),bybd(2)/sqrt(rhobd(2)),
     +            sqrt(0.5*gamma*pbd(2)/rhobd(2))
c
      return
      end
