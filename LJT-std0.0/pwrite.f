	subroutine pwrite
c***********************************************************************
	include 'misflin'
c
c       Version 2p
c       This version includes code which generates output pertaining
c       to the `ballistic relaxation.'
c
	integer  i, j, ix, iy, iz
c.....................................................
      write(26,10)  istep
      write(26,11)  nx,ny,nz,
     +               intart,isafe,iend,nsmooth,ivisrho,ivisu,
     +               igrid,idiag,iferror,
     +               start,ieein,ieeout,jout,fluxcal,mitsat,deltj
      write(26,111)  relax,nrelax,nmrelax

      write(26,12)  dt,gamma,eta,visx,visy,visz,
     +               xmin,ymin,zmin,xmax,ymax,zmax
      if (etasw.eq.0) write(26,120)
      if (etasw.eq.1) write(26,121)
      if (etasw.ne.0 .and. etasw.ne.1) write(26,122) etasw,jcrit0
      write(26,13)  psi,phi,bmsh,pmsp,kappa,
     +              rho0,delrho,xrho,dxrho,by0,bz0,delbz,
     +              bzmsp,bzmsh,pmsp,pmsh,betasp,betash,
     +              vamsp,vamsh,csmsp,csmsh
     
      write(26,131)  nsim,bsim,vsim,psim,lsim,tsim
      write(26,132)  rhobd(1),pbd(1),vxbd(1),vybd(1),vzbd(1),
     +                   bxbd(1),bybd(1),bzbd(1),
     +                   rhobd(2),pbd(2),vxbd(2),vybd(2),vzbd(2),
     +                   bxbd(2),bybd(2),bzbd(2)
     
      write(26,14)  aequi(1),aequi(2),aequi(3),
     +               zentr(1),zentr(2),zentr(3),
     +               eps(1),eps(2),eps(3)
      write(26,15)  perio(1),perio(2),perio(3),
     +               lsym(1,1),lsym(2,1),lsym(1,2),lsym(2,2),
     +               lsym(1,3),lsym(2,3)
      do 100 j = 1,3
      do 100 i = 1,2
	 write(26,16) i,j,
     +                 kbx(i,j),kby(i,j),kbz(i,j),
     +                 ksx(i,j),ksy(i,j),ksz(i,j),
     +                 arho(i,j),au(i,j),ares(i,j),
     +                 abx(i,j),aby(i,j),abz(i,j),
     +                 asx(i,j),asy(i,j),asz(i,j)
  100 continue

      do 101 j = 1,3
      do 101 i = 1,2
	write(26,17) i,j,
     +    crho1(i,j),crho2(i,j),cu1(i,j),cu2(i,j),
     +    cres1(i,j),cres2(i,j),
     +    cbx1(i,j),cbx2(i,j),cby1(i,j),cby2(i,j),cbz1(i,j),cbz2(i,j),
     +    csx1(i,j),csx2(i,j),csy1(i,j),csy2(i,j),csz1(i,j),csz2(i,j)
  101 continue
c
 1003   continue
        write(26,8891)
        do 1011 ix = 1,nx
 	 write(26,8881) x(ix),difx(ix),ddifx(ix),ddifpx(ix),ddifmx(ix),
     +                  meanpx(ix),meanmx(ix)
 1011   continue
        write(26,8892)
        do 1012 iy = 1,ny
 	 write(26,8881) y(iy),dify(iy),ddify(iy),ddifpy(iy),ddifmy(iy),
     +                  meanpy(iy),meanmy(iy)
 1012   continue
        write(26,8893)
        do 1013 iz = 1,nz
	 write(26,8881) z(iz),difz(iz),ddifz(iz),ddifpz(iz),ddifmz(iz),
     +                  meanpz(iz),meanmz(iz)
 1013   continue
c
   10 format(/1x,'program parameter(intstep=',i5,') : ',
     +       /1x,'---------------------------------------')
c
   11 format(/1x,'nx     =',i5,3x,'ny      =',i5,3x,'nz      =',i5,3x,
     +       /1x,'intart =',i5,3x,'isafe   =',i5,3x,'iend    =',i5,3x,
     +       /1x,'nsmooth=',i5,3x,'ivisrho =',i5,3x,'ivisu   =',i5,3x,
     +       /1x,'igrid  =',i5,3x,'idiag   =',i5,3x,'iferror =',i5,3x,
     +       /1x,'start  =',l5,3x,'ieein   =',l5,3x,'ieeout  =',l5,3x,
     +       /1x,'jout   =',l5,3x,'fluxcal =',l5,3x,'mitsat  =',l5,3x,
     +       /1x,'deltj  =',f8.2)
c
  111 format(/1x,'Relaxation:',l5,3x,
     +           'Every',i5,2x,'time steps up to',i7)
c
   12 format(/1x,'timestep =',f9.5,3x,'gamma =',f9.5,
     +        3x,'resistivity =',f9.5,
     +       /1x,'visx =',f9.5,3x,'visy =',f9.5,3x,'visz =',f9.5,
     +       /1x,'xmin =',f9.3,3x,'ymin =',f9.3,3x,'zmin =',f9.3,
     +       /1x,'xmax =',f9.3,3x,'ymax =',f9.3,3x,'zmax =',f9.3)
  120 format(/1x,'resistivity is constant')
  121 format(/1x,'resistivity is time dependent')
  122 format(/1x,'resistivity is parameter dependent; model:',i2,
     +           '   crit current density:', f8.3)
c
c
   13 format(/1x,' equilibrium parameter : ',
     +       /1x,'---------------------------',
     +       /1x,'phi    =',f7.2,3x,'psi   =',f7.2,3x,
     +       /1x,'bmsh   =',f7.2,3x,'pmsp  =',f7.2,3x,'kappa  =',f7.2,
     +       /1x,'rho0   =',f7.2,3x,'delrho=',f7.2,3x,
     +           'xrho   =',f7.2,3x,'dxrho  =',f7.2,
     +       /1x,'which yield:  by0 =',f8.3,3x,'  bz0 =',f8.3,3x,
     +           'delbz =',f8.3,
     +       /1x,'                !      msp      !      msh      ! ',
     +       /1x,'--------------------------------------------------',
     +       /1x,'  bz            ! ',f12.5,'  ! ',f12.5,'  !',
     +       /1x,'  p             ! ',f12.5,'  ! ',f12.5,'  !',
     +       /1x,'  beta          ! ',f12.5,'  ! ',f12.5,'  !',
     +       /1x,'  v-alfven      ! ',f12.5,'  ! ',f12.5,'  !',
     +       /1x,'  v-schall      ! ',f12.5,'  ! ',f12.5,'  !'  )
     
  131 format(/1x,'Normalisation for : ',
     +       /1x,' No density ! magn. field!  velocity  !  pressure  !',
     +           'length units!    time     ',
     +       /1x,'  cm**(-3)  !     nT     !    km/s    !   nPascal  !',
     +           '     km     !      s      ',
     +         /,f11.2,f13.2,f13.2,f13.4,f13.2,f13.4)
  132 format(/1x,'Initial asymptotic values at xmin (1. row) ',
     +           'and xmax (2. row) for: ',
     +       /1x,'    rho      p        vx       vy       vz    ',
     +           '   bx       by       bz    ',
     +       /,8f9.3/,8f9.3)
     
   14 format(/1x,'uniform grid  in       x:',l3,8x,'y:',l3,
     +        8x,'z:',l3,
     +       /1x,'max resolution centered for x:',l3,
     +        8x,'y:',l3,8x,'z:',l3,
     +       /1x,'with grid spacing:         x:',f7.3,
     +        4x,'y:',f7.3,4x,'z:',f7.3)

   15 format(/1x,'periodic boundary conditions in     x:',l4,3x,
     +        'in y:',l4,3x,'in z:',l4,
     +       /1x,'line symmetre (entlang z) at xmin:',l4,
     +        3x,'xmax:',l4,
     +       /1x,'                (entlang x) at ymin:',l4,
     +        3x,'ymax:',l4,
     +       /1x,'                (entlang x) at zmin:',l4,
     +        3x,'zmax:',l4)

   16 format(/1x,'min/max = ',i4,5x,'x/y/z = ',i4,
     +       /1x,'kbx  =',i4,3x,'kby  =',i4,3x,'kbz  =',i4,3x,
     +       /1x,'ksx  =',i4,3x,'ksy  =',i4,3x,'ksz  =',i4,3x,
     +       /1x,'arho =',f5.2,2x,'au   =',f5.2,2x,'ares =',f5.2,2x,
     +       /1x,'abx  =',f5.2,2x,'aby  =',f5.2,2x,'abz  =',f5.2,2x,
     +       /1x,'asx  =',f5.2,2x,'asy  =',f5.2,2x,'asz  =',f5.2,2x)
c
   17 format(/1x,'min/max = ',i4,5x,'x/y/z = ',i4,
     +       /1x,'crho1=',f5.2,2x,'crho2=',f5.2,2x,
     +           'cu1  =',f5.2,2x,'cu2  =',f5.2,2x,
     +           'cres1=',f5.2,2x,'cres2=',f5.2,2x,
     +       /1x,'cbx1 =',f5.2,2x,'cbx2 =',f5.2,2x,
     +           'cby1 =',f5.2,2x,'cby2 =',f5.2,2x,
     +           'cbz1 =',f5.2,2x,'cbz2 =',f5.2,2x,
     +       /1x,'csx1 =',f5.2,2x,'csx2 =',f5.2,2x,
     +           'csy1 =',f5.2,2x,'csy2 =',f5.2,2x,
     +           'csz1 =',f5.2,2x,'csz2 =',f5.2)
c
 8881 format(1x,7f10.3)
 8891 format(/4x,'  x         dif       ddif      ddifp     ddifm',
     +           '     meanp     meanm ')
 8892 format(/4x,'  y         dif       ddif      ddifp     ddifm',
     +           '     meanp     meanm ')
 8893 format(/4x,'  z         dif       ddif      ddifp     ddifm',
     +           '     meanp     meanm ')
c
      return
      end
