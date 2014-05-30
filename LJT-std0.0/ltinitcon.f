	subroutine initcon
c***********************************************************************
	include 'misflin'
c
c-------for array
	integer  ix,iy,iz
c-------for MHD 
	real     gaminv,eta0,ptot,mu0
        real     localx(nx)
c------------for Jupiter geometry
        real    aa,bb,cc,xj,rr

       

c ...................................................
c  initial configuration
c ...................................................
      eta0   = 0.0
      help   = 0.0
 
      if (etasw.eq.0)  eta0 = eta
      gaminv = 1./gamma



c------------geometry for Jupiter
c--   ((x-xj)/aa)^2+((y-0.0)/bb)^2+((z-0.)/c)^2=R^2



            xj=15 
            aa = 1.0
            bb = 2.0
            cc = 1.0
            rr = 1.5
            rr=rr*rr
            


c-----Solar wind in x dirtion 
            
            ptot=pmsp+by0*by0+bz0*bz0
            
           
           
        help(:,1,1) = 0.5*(1-tanh(x-0.75*(xj-xmin)))  
        rhoprof = rho0+delrho*help(:,1,1)
         bxprof = 0.
         byprof = by0*help(:,1,1)
         bzprof = bz0*help(:,1,1)

         uprof = ptot-bzprof*bzprof-byprof*byprof
	 uprof = (0.5*uprof)**gaminv

          sxprof = vx0*help(:,1,1)*rhoprof
	  syprof = 0.
	  szprof = 0. 
          help   = 0.
         
        
 
               


c------------

      
      if (start .eqv. .false.) return 

         
  

c-----------------------3d set


            mu0 = 10.


        do 60 iz = 1,nz
	do 60 iy = 1,ny
	do 60 ix = 1,nx
 	   rho(ix,iy,iz) = rhoprof(ix)
	     u(ix,iy,iz) =   uprof(ix)
            bx(ix,iy,iz)  = bxprof(ix)
            by(ix,iy,iz)  = byprof(ix)
            bz(ix,iy,iz)  = bzprof(ix)
    	    sx(ix,iy,iz)  = sxprof(ix)
   	    sy(ix,iy,iz)  = syprof(ix)
            sz(ix,iy,iz)  = szprof(ix)
          prof(ix,iy,iz) = 1./cosh(x(ix))/cosh(y(iy))/cosh(z(iz))
           res(ix,iy,iz) = eta0*prof(ix,iy,iz)

c---------------
        muprof(ix,iy,iz)= ((x(ix)-xj)/aa)**2+(y(iy)/bb)**2+(z(iz)/cc)**2
        muprof(ix,iy,iz)= mu0/cosh(muprof(ix,iy,iz)/rr)
c        muprof(ix,iy,iz)= mu0*(1-tanh(muprof(ix,iy,iz)-rr))

   60 continue

      call binout(0)


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





