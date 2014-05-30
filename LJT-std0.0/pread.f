	subroutine pread
c***********************************************************************
	include 'misflin'
c
c       Version 2p
c       This version includes the integer variables nrelax and nmrelax
c       among the variables read from the file magin.  These variables
c       are used for the `ballistic relaxation.'
c
	logical  ldum
	integer  idum, i, j
	real     rdum
c.....................................................
      read (15,2) etasw,isafe,iend,intart,ivisrho,ivisu,nsmooth,idiag,
     +            nrelax,nmrelax
      write (6,2) etasw,isafe,iend,intart,ivisrho,ivisu,nsmooth,idiag,
     +            nrelax,nmrelax
      if (start) then
        read (15,2) nzz
        read (15,3) cval
      else
        read (15,2) idum
        read (15,3) rdum
      end if

      read (15,3) deltj,dt,gamma,eta,jcrit0,visx,visy,visz
      write (6,3) deltj,dt,gamma,eta,jcrit0,visx,visy,visz

      if (start) then
	read (15,4) xmin,xmax,ymin,ymax,zmin,zmax
        read (15,4) psi,phi,bmsh,pmsp,kappa,rho0,delrho,xrho,dxrho
        read (15,4) p0,bx0,by0,bz0,vx0,vy0,vz0
        read (15,4) rho1,p1,bx1,by1,bz1,vx1,vy1,vz1
        read (15,4) nsim,bsim,lsim
	read (15,1) aequi(1),aequi(2),aequi(3)
        read (15,1) zentr(1),zentr(2),zentr(3)
        read (15,3) eps(1),eps(2),eps(3)
        write (6,4)  xmin,xmax,ymin,ymax,zmin,zmax
        write (6,4)  phi,bmsh,pmsp,kappa,rho0,delrho,xrho,dxrho
        write (6,1)  aequi(1),aequi(2),aequi(3)
        write (6,1)  zentr(1),zentr(2),zentr(3)
        write (6,3)  eps(1),eps(2),eps(3)
      else
	read (15,4) rdum,rdum,rdum,rdum,rdum,rdum
	read (15,4) rdum,rdum,rdum,rdum,rdum,rdum,rdum,rdum,rdum
	read (15,4) rdum,rdum,rdum,rdum,rdum,rdum,rdum
	read (15,4) rdum,rdum,rdum,rdum,rdum,rdum,rdum,rdum
	read (15,4) rdum,rdum,rdum
	read (15,1) ldum,ldum,ldum
	read (15,1) ldum,ldum,ldum
	read (15,3) rdum,rdum,rdum
      end if

      if (start)  then
        read (15,1)  perio(1),perio(2),perio(3)
        write (6,1)  perio(1),perio(2),perio(3)
        read (15,1)  lsym(1,1),lsym(2,1),lsym(1,2),lsym(2,2),
     +               lsym(1,3),lsym(2,3)
        write (6,1)  lsym(1,1),lsym(2,1),lsym(1,2),lsym(2,2),
     +               lsym(1,3),lsym(2,3)
        do 50 j = 1,3
        do 50 i = 1,2
	   read(15,2) kbx(i,j),kby(i,j),kbz(i,j),
     +                 ksx(i,j),ksy(i,j),ksz(i,j)
	   read(15,3) arho(i,j),au(i,j),ares(i,j),
     +                 abx(i,j),aby(i,j),abz(i,j),
     +                 asx(i,j),asy(i,j),asz(i,j)
   50   continue
      else
	read (15,1) ldum,ldum,ldum
	read (15,1) ldum,ldum,ldum,ldum,ldum,ldum
        do 55 j = 1,6
          read (15,2) idum,idum,idum,idum,idum,idum
	  read (15,3) rdum,rdum,rdum,rdum,rdum,rdum,rdum,rdum,rdum
   55   continue
      end if
      
      if (mitsat) then
        read(15,4) tsat0,tsatout 
	do 100 i = 1,npos
          read(15,4) xpos(i),ypos(i),zpos(i),  
     +               dxpos(i),dypos(i),dzpos(i),
     +               vxpos(i),vypos(i),vzpos(i)
  100   continue
      end if
      
c   normaliszation of simulation, nsim, bsim, lsim from parlist
c                              in cm**(-3), nT, and km
c   vsim, psim, tsim in km/s, nPascal, and s
      vsim = 21.8*bsim/sqrt(nsim)
      psim = 0.01*bsim**2/8.0/pi
      tsim = lsim/vsim

c
    1 format(1x,l10)
    2 format(1x,i9)
    3 format(1x,f14.8)
    4 format(1x,f14.5)
c
      return
      end
