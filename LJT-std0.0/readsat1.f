	program   readsat
c
        integer    nsat,nt,nc
        real       pi
	parameter  (pi = 3.14159265536)
	parameter  (nsat = 60, nt = 3000, nc = 4)
c
	integer  i,k,ic(nc),ntmax
	real     zeitsat,xsat(nsat),ysat(nsat),zsat(nsat),
     +           vxsat(nsat),vysat(nsat),vzsat(nsat),
     +           bxs(nsat),bys(nsat),bzs(nsat),
     +           vxs(nsat),vys(nsat),vzs(nsat),
     +           rhosat(nsat),psat(nsat),bsat(nsat),
     +           nnorm,bnorm,vnorm,pnorm,lnorm,tnorm,
     +           rhobd(2),pbd(2),vxbd(2),vybd(2),vzbd(2),
     +           bxbd(2),bybd(2),bzbd(2),
     +           yd(nc,nt),zd(nc,nt),
     +           byd(nc,nt),bzd(nc,nt),vyd(nc,nt),vzd(nc,nt),
     +           xc(nc,nt),yc(nc,nt),zc(nc,nt),
     +           bxc(nc,nt),byc(nc,nt),bzc(nc,nt),
     +           vxc(nc,nt),vyc(nc,nt),vzc(nc,nt),
     +           rhoc(nc,nt),pc(nc,nt),bc(nc,nt),ptotc(nc,nt),
     +           time,zeit(nt),phi,cphi,sphi,
     +           rhom(2),pm(2),vxm(2),vym(2),vzm(2),
     +           bxm(2),bym(2),bzm(2),
     +           xsatm(nsat),ysatm(nsat),zsatm(nsat),
     +           vxsatm(nsat),vysatm(nsat),vzsatm(nsat)
c This programm reads data for numerical probes in 2 and 3D simulation.
c For nc(=4) it performs a coordinate + Galilei transformation into the 
c magnetospheric frame where the magnetospheric field is along z, 
c x is normal to the initial current layer (pointing into the magnetosheath, 
c and the plasma is at rest. 
c It also applies a reasonable normalization. 

c
c---READ SATALLITE PARAMETERS FROM FILE SATBIN 
c 
c    !!!ALL DATA IS HERE IN SIMULATION COORDINATES AND NORMALIZED!!!
 
      open (39,file='satbin',form='unformatted')
c   For Cray comment above line and uncomment following line
c      call asnunit(39,'-a satbin -F f77 -N ieee',ier)
c 1. nssat: number of probes/datasets=nsat!!!
c    zeitsat: startime when first data is recorded
      read(39) nssat,zeitsat
      if ( nssat.ne.nsat) then
	  write(31,9)  nsat, nssat
	  write(*,9)  nsat, nssat
	  stop
      end if
c 2. next read inputs suggested normalization parameters

      read(39) nnorm,bnorm,vnorm,pnorm,lnorm,tnorm
      
c 3. Asymptotic density,pressure, velocity, and magn. field 
c    at the boundaries xmin (MSP, index 1) and xmax (MSH, index 2)

      read(39) rhobd,pbd,vxbd,vybd,vzbd,bxbd,bybd,bzbd
      
c 4. Initial locations and velocities of probes; velocity is relative to 
c    plasma velocity at xmin, i.e., measured in the restframe of
c    the magnetospheric boundary 

      read(39) xsat,ysat,zsat,vxsat,vysat,vzsat
c
c---To change NORMALIZATION uncomment the following lines
c     Normaliszation: nnorm(cm**(-3)), bnorm(nT), lnorm(km) 
c                 =>  vnorm(km/s), pnorm(nPascal), tnorm(s) 
c   E.g.:
c     nnorm = 20.0
c     bnorm = 10.0
c     lnorm = 200.0
c     vnorm = 21.8*bnorm/sqrt(nnorm)
c     pnorm = 0.01*bnorm**2/8.0/pi
c     tnorm = lnorm/vnorm
c
c---WRITE SATALLITE PARAMETERS TO ASCII FILE(SATDAT1) AND CONSOLE   
      open (31,file='satdat1')
      write(31,21)  nssat,zeitsat
      write(*,21)   nssat,zeitsat
      write(31,210) 
      write(*,210) 
      write(31,211) nnorm,bnorm,vnorm,pnorm,lnorm,tnorm
      write(*,211)  nnorm,bnorm,vnorm,pnorm,lnorm,tnorm
      write(31,212)  rhobd(1),pbd(1),vxbd(1),vybd(1),vzbd(1),
     +                bxbd(1),bybd(1),bzbd(1),
     +                rhobd(2),pbd(2),vxbd(2),vybd(2),vzbd(2),
     +                bxbd(2),bybd(2),bzbd(2)
      write(*,212)  rhobd(1),pbd(1),vxbd(1),vybd(1),vzbd(1),
     +                bxbd(1),bybd(1),bzbd(1),
     +                rhobd(2),pbd(2),vxbd(2),vybd(2),vzbd(2),
     +                bxbd(2),bybd(2),bzbd(2)
      write(31,22) 
      write(*,22) 
      do 60 i = 1,nsat
          write(31,23) i,xsat(i),ysat(i),zsat(i),
     +                   vxsat(i),vysat(i),vzsat(i)
          write(*,23) i,xsat(i),ysat(i),zsat(i),
     +                   vxsat(i),vysat(i),vzsat(i)
   60 continue
      write(31,24) 
      write(*,24) 
c---SELECT NC=4 SATALLITE INDICES  
      write(*,*) 'Choose 4 satellite indices (<=',nsat,'). 1. index?'
      read(*,*) ic(1)
      write(31,241) ic(1)
      do 65 i = 2, nc
        write(*,*) i,'. index?'
        read(*,*) ic(i)
        write(31,242) ic(i)
   65 continue
       
c---READ DATA AND TRAJECTORIES FOR SELECTED INDICES   
      do 70 i = 1,nt
       ntmax = i-1
       read(39,end=100) time,bxs,bys,bzs,vxs,vys,vzs,rhosat,psat,bsat
c     +                  ,xsat,ysat,zsat
c   NOTE: here xsat,ysat,zsat are the satellite coordinates in the 
c         simulation box which have an additional velocity of 
c         vxbd(1),vybd(1),vzbd(1) (i.e.at the magnetospheric boundary) 
c         in order to be approximately at rest in the frame of the 
c         magnetosphere. The simulation box is moving with half 
c         the magnetosheath plasma velocity for Kelvin Helmholtz cases!
c         All following parameters with 'label' c are in the restframe 
c         of the magnetosphere (defined by initial plasma flow at the 
c         magnetospheric boundary), but still in simulation coordinates
c
       zeit(i) = tnorm*time
       do 70 k = 1,nc
         xc(k,i)  = lnorm*( xsat(ic(k)) + (time-zeitsat)*vxsat(ic(k)) )
         yd(k,i)  = lnorm*( ysat(ic(k)) + (time-zeitsat)*vysat(ic(k)) )
         zd(k,i)  = lnorm*( zsat(ic(k)) + (time-zeitsat)*vzsat(ic(k)) )
         bxc(k,i) = bnorm*bxs(ic(k))
         byd(k,i) = bnorm*bys(ic(k))
         bzd(k,i) = bnorm*bzs(ic(k))
         vxc(k,i) = vnorm*vxs(ic(k))
         vyd(k,i) = vnorm*vys(ic(k))
         vzd(k,i) = vnorm*vzs(ic(k))
         rhoc(k,i)= nnorm*rhosat(ic(k))
         pc(k,i)  = pnorm*psat(ic(k))
         bc(k,i)  = bnorm*bsat(ic(k))
         ptotc(k,i)= pnorm*( psat(ic(k)) + bsat(ic(k))**2 )
   70 continue
  100 write(*,*) 'last record:', ntmax
      close (39)
      write(31,243) ntmax
c---Rotate y and z compents to the assumed magnetospheric coordinate system
c     with z locally northward
c
c - first determine rotation angle:
c      
       write(*,*) '  INPUT Rotation Angle, PHI.'
       write(*,*) '  For KH type a negative value'
       write(*,*) '              -> PHI  = atan(bybd(1)/bzbd(1)).'
       write(*,*) '  For 2D Reconnection and pressure pulse cases'
       write(*,*) '              PHI ought to be typically 90 degrees.'
       read(*,*) phi
       phi = pi*phi/180.
       if (phi .lt. 0) phi  = atan(bybd(1)/bzbd(1))
       sphi = sin(phi)
       cphi = cos(phi)
       write(31,325) 180*phi/pi
       write(*,325) 180*phi/pi
       write(*,*) '  All output is recorded in the file SATDAT1'
c - then satellite data and trajectories:
       do 80 k = 1,nc
       do 80 i = 1,ntmax
         yc(k,i)  = yd(k,i)*cphi - zd(k,i)*sphi
         zc(k,i)  = yd(k,i)*sphi + zd(k,i)*cphi
         byc(k,i) = byd(k,i)*cphi - bzd(k,i)*sphi
         bzc(k,i) = byd(k,i)*sphi + bzd(k,i)*cphi
         vyc(k,i) = vyd(k,i)*cphi - vzd(k,i)*sphi
         vzc(k,i) = vyd(k,i)*sphi + vzd(k,i)*cphi
   80 continue
c - initial asymptotic values in the magnetospheric frame  
c   with reasonable units
      do 85 i = 1,2
        rhom(i) = nnorm*rhobd(i)
        pm(i) = pnorm*pbd(i)
        vxm(i) = vnorm*vxbd(i)
        vym(i) = vnorm*( vybd(i)*cphi - vzbd(i)*sphi )
        vzm(i) = vnorm*( vybd(i)*sphi + vzbd(i)*cphi )
        bxm(i) = bnorm*bxbd(i)
        bym(i) = bnorm*( bybd(i)*cphi - bzbd(i)*sphi )
        bzm(i) = bnorm*( bybd(i)*sphi + bzbd(i)*cphi )
   85 continue
      vym(2) = vym(2) - vym(1)
      vym(1) = 0
      vzm(2) = vzm(2) - vzm(1)
      vzm(1) = 0
c - sat locations and velocity in the magnetospheric frame  
      do 90 i = 1,nsat
        xsatm(i) = lnorm*xsat(i)
        ysatm(i) = lnorm*( ysat(i)*cphi - zsat(i)*sphi )
        zsatm(i) = lnorm*( ysat(i)*sphi + zsat(i)*cphi )
        vxsatm(i) = vnorm*vxsat(i)
        vysatm(i) = vnorm*( vysat(i)*cphi - vzsat(i)*sphi )
        vzsatm(i) = vnorm*( vysat(i)*sphi + vzsat(i)*cphi )
   90 continue
      write(31,326)  rhom(1),pm(1),vxm(1),vym(1),vzm(1),
     +                bxm(1),bym(1),bzm(1),
     +                rhom(2),pm(2),vxm(2),vym(2),vzm(2),
     +                bxm(2),bym(2),bzm(2)
      write(*,326)  rhom(1),pm(1),vxm(1),vym(1),vzm(1),
     +                bxm(1),bym(1),bzm(1),
     +                rhom(2),pm(2),vxm(2),vym(2),vzm(2),
     +                bxm(2),bym(2),bzm(2)
      write(31,327) 
      write(*,327) 
      do 95 i = 1,nsat
          write(31,328) i,xsatm(i),ysatm(i),zsatm(i),
     +                   vxsatm(i),vysatm(i),vzsatm(i)
          write(*,328) i,xsatm(i),ysatm(i),zsatm(i),
     +                   vxsatm(i),vysatm(i),vzsatm(i)
   95 continue

c---WRITE DATA AND TRAJECTORIES TO ASCII FILE SATDAT1
      do 170 k = 1,nc
       write(31,251) ic(k)
       write(31,25) 
       do 170 i = 1,ntmax
	 write(31,26) zeit(i),bxc(k,i),byc(k,i),bzc(k,i),
     +                        vxc(k,i),vyc(k,i),vzc(k,i),
     +                        rhoc(k,i),pc(k,i),bc(k,i),ptotc(k,i)
  170 continue
      do 180 k = 1,nc
       write(31,27) ic(k)
       write(31,271) 
       do 180 i = 1,ntmax
	 write(31,28) zeit(i), xc(k,i), yc(k,i), zc(k,i)
  180 continue

      close (31)

    9 format(1x,'!!!! Program parameter nsat(',i4,') is wrong !!!!!',
     +       /1x,'!!!!    Change parameter to:',i4,'    !!!!!')
   21 format(1x,'Number of satellites       start time',/1x,i10,f22.3)
  210 format(/1x,'SIMULATION PARAMETERS')
  211 format(/1x,'Normalisation for : ',
     +       /2x,' No density ! magn. field!  velocity  !  pressure  !',
     +           'length units!    time     ',
     +       /2x,'  cm**(-3)  !     nT     !    km/s    !   nPascal  !',
     +           '     km     !      s      ',
     +         /,f11.2,f13.2,f13.2,f13.4,f13.2,f13.4)
  212 format(/1x,'Initial asymptotic values at xmin',
     +           ' (1. row, magnetosphere) ',
     +       /1x,' and xmax (2. row) in the simulation frame: ',
     +       /1x,'    rho      p        vx       vy       vz    ',
     +           '   bx       by       bz    ',
     +       /,8f9.3/,8f9.3)
   22 format(/1x,'Initial satellite locations and velocity:',
     +       /1x,'sat  !    x    !    y    !    z    !',
     +                 '    vx   !    vy   !    vz   ')
   23 format(1x,i3,'   ',f7.2,'   ',f7.2,'   ',f7.2,
     +             '   ',f7.2,'   ',f7.2,'   ',f7.2)
   24 format(1x,'Satellite velocity is relative to restframe of the',
     +          ' magnetosphere (at xmin),',
     +       /1x,' however, still in simulation coordinates!!!!!!')
  325 format(/1x,'PARAMETERS IN MAGNETOSPHERIC FRAME and COORDINATES',
     +       /1x,'  Magnetospheric Coordinates defined such that ',
     +           'the magnetic field at the ',
     +       /1x,'  magnetospheric boundary is in the ',
     +            'positive z direction',
     +       /1x,'  This requires a rotation of y and z axes by',
     +            f7.2,' degrees')
  326 format(/1x,'Initial asymptotic values at xmin (1. row)',
     +       /1x,' and xmax (2. row) in the magnetospheric frame: ',
     +       /1x,'    rho      p        vx       vy       vz    ',
     +           '   bx       by       bz    ',
     +       /1x,' cm**(-3)   nP       km/s     km/s     km/s   ',
     +           '   nT       nT       nT    ',
     +       /,f9.2,f9.3,6f9.2,
     +       /,f9.2,f9.3,6f9.2)
  327 format(/1x,'Initial satellite locations and velocity in magnet. ',
     +           ' frame:',
     +       /1x,'sat  !    x     !    y     !    z     !',
     +                 '    vx    !    vy    !    vz    ',
     +       /1x,'     !    km    !    km    !    km    !',
     +                 '   km/s   !   km/s   !   km/s   ')
  328 format(1x,i3,'   ',f8.1,'   ',f8.1,'   ',f8.1,
     +             '   ',f8.1,'   ',f8.1,'   ',f8.1)
  241 format(/1x,'Chosen satellite indices: ',i4)
  242 format(1x,i30)
  243 format(1x,'Size of datasets:', i4)
   25 format(1x,' time    bx     by     bz     vx     vy     vz   ',
     +          'rhosat  psat    bsat   ptot',
     +      /1x,'   s     nT     nT     nT    km/s   km/s   km/s  ',
     +          'cm**-3   nP     nT     nP ')
  251 format(/1x,'Satellite data for Index:',i4)
   26 format(1x,f6.2,3f7.2,3f7.1,f7.2,f7.3,f7.2,f7.3)
   27 format(/1x,'Satellite trajectory for Index:',i4)
  271 format(1x,'  time        x           y           z')
   28 format(1x,f6.2,3f12.3)

      stop
      end
