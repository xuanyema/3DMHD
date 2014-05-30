	program   magnfl
c***********************************************************************
	include 'misflin'
c
	integer          ksafe,kdia,msafe,mout
c................................................................
c  Note: 1. To run program with output of current density or 
c           flux surface set parameters in magin, modify makefile to 
c           compile additional subroutines, and uncomment lines following 
c           markers  current! or flux! in: magnfl.f, orgstart.f, misflin
c                                          
c        2. To run program on PC or workstation set ieein/ieeout to false 
c           and comment lines containing asnunit or assign in: 
c           binin.f, binout.f, and satstart.f
c................................................................
c   i) bookkeeping of program, open und close, formatted output
c................................................................
      open (15,file='magin')
      open (17,file='magdia1')
      open (18,file='magdia2')
      open (19,file='magdfte')
      open (25,file='magdmax')
      open (26,file='magpar')
      open (27,file='magdfl1')
      open (28,file='magdfl2')
c
c   files im programm:
c     14=magtapn, 15=magin,   
c     17=magdia1, 18=magdia2, 19=magdfte, 
c     25=magdmax, 26=magpar, 
c     27=magdfl1, 28=magdfl2, 
c     31=sat1,    32=sat2,    33=sat3, 
c     34=sat4,    35=sat5,    36=sat6,    39=satbin
c
      call orgstart
      if (fini) stop
      msafe=0
      mout=0
      ldiag = ldiag + 1
      tdiag(ldiag) = time
      call diag0
      call diag1
      call diag2
c uncomment the following lines for 
c   flux! computation or
c   current! density output
c      if (time .ge. 50.) then
c       if ( fluxcal )  call orgflux
c       if ( jout )  call orgjout
c      endif

c    -----------------------------------
c   ii)  start of integration (first 2 steps
c    -----------------------------------
  400 call intstart
      if (nzz .lt. nz) call nzznew
      if ( ferror ) go to 900
      if ( istep .ge. iend) goto 900
c    -----------------------------------
c     iii) diagnostics and output
c    -----------------------------------
  500 if (mod(istep,100).eq.0) write(26,58) istep
c    -----------------------------------
c     a) diagnostics 
c    -----------------------------------
      kdia = mod(istep,idiag)
      if ( kdia .eq. 0 )  then
	if (ldiag .ge. ndiag) then
	    call diagout
	    ldiag = 0
	end if
	ldiag = ldiag + 1
	tdiag(ldiag) = time
	call diag0
	call diag1
	call diag2
      end if
      if (mod(istep,10*idiag) .eq. 0) call diag3
c    -----------------------------------
c   flux integration and current density output
c    -----------------------------------
c uncomment the following lines for 
c   flux! computation or
c   current! density output
c      if (time .ge. 50.) then
c        if ( fluxcal .and. (tflux .ge. deltfl) )  then
c         call orgflux
c         tflux = tflux-deltfl
c        endif
c        tflux = tflux+2.*dt
c        if ( jout .and. (tjout .ge. deltj) )  then
c         call orgjout
c         tjout = tjout-deltj
c        endif
c        tjout = tjout+2.*dt
c      endif
c    -----------------------------------
c     b) satellites
c    -----------------------------------
      if ( time .ge. tsat0 .and. mitsat)  then
        if (.not. startsat) call satstart
        call satdat
      endif
c    -----------------------------------
c     c) binary output (ksafe)
c    -----------------------------------
      ksafe = mod(istep,isafe)
      if ( ksafe.eq.0 .or. newdt) then
         igrid=1
	 call lax
	 call bound
	 call termin
	 if(newdt) dt = 0.9*dt
	 newdt = .false.
	 if ( ferror .and. iferror.le.nsmooth ) call osmooth
	 if ( ferror ) go to 900
	 if ( ksafe .eq. 0 ) then
	   msafe=msafe+1
	   call binout(msafe)
	   write(26,52) istep
	 end if
      end if
      if ( istep .ge. iend) goto 900
      if ( ksafe.eq.0 ) goto 400
c
      call orgint
      if (nzz .lt. nz) call nzznew
      if ( ferror ) goto 900
c
      goto 500
  900 call diagout
c

      close (15)
      close (17)
      close (18)
      close (19)
      close (25)
      close (26)
      close (27)
      close (28)
      if ( startsat)  then
        call assign('assign -R')
        close (49)
        close (31)
        close (32)
        close (33)
        close (34)
        close (35)
        close (36)
      endif
    1 format(1x,l10)
   52 format(1x,'binary data written',i5)
   58 format(1x,'at integration step:',i6)
c
c  list of subroutines:
c ------------------------
c main magn:
c block data vorbes0:
c        orgstart:
c        intstart:
c        orgint:
c        oglatt:
c        lespar:
c        schreib:
c        grid:
c        anfang:
c        orgaus:
c        aus:
c        glatt:
c        meanlax:
c        lax:
c        leap:
c        rdcoef:
c        rand:
c        symfx:
c        symafx:
c        symtx:
c        symty:
c        symtz:
c        abbruch:
c        diag1:
c        diag2:
c        diagaus:
c
      stop
      end
	block data first
c***********************************************************************
	include 'misflin'
c.................................................................
      data    mass/ndiag*0./,massu/ndiag*0./,
     +        enmag/ndiag*0./,enkin/ndiag*0./, enthe/ndiag*0./,
     +        vgradp/ndiag*0./, vdjxb/ndiag*0./, resj2/ndiag*0./,
     +        ekinpu/ndiag*0./, ethpu/ndiag*0./, ebpu/ndiag*0./,
     +        frho/ndiag*0./, fekin/ndiag*0./, fethe/ndiag*0./,
     +        feb/ndiag*0./, febi/ndiag*0./, febr/ndiag*0./,
     +        fphi/ndiag*0./, fphip/ndiag*0./, fphipi/ndiag*0./,
     +        fphipr/ndiag*0./
      data    isat/40/,timesat/0.0/,
     +        xpos/npos*0./,ypos/npos*0./,zpos/npos*0./,
     +        dxpos/npos*0./,dypos/npos*0./,dzpos/npos*0./,
     +        vxpos/npos*0./,vypos/npos*0./,vzpos/npos*0./,
     +        xsat/nsat*0./,ysat/nsat*0./,zsat/nsat*0./,
     +        vxsat/nsat*0./,vysat/nsat*0./,vzsat/nsat*0./,
     +        bxs/nsat*0./,bys/nsat*0./,bzs/nsat*0./,
     +        vxs/nsat*0./,vys/nsat*0./,vzs/nsat*0./,
     +        rhosat/nsat*0./,psat/nsat*0./,bsat/nsat*0./,
     +        ptotsat/nsat*0./
      data    fini/.false./,start/.true./,etasw/0/,
     +        mitsat/.true./,ferror/.false./,newdt/.false./,
     +        istep/0/,ieein/.true./,isafe/800/,iend/8000/,intart/1/,
     +        ivisrho/1200/,ivisu/40000/,nsmooth/80/,ieeout/.true./,
     +        iferror/0/,igrid/1/,
     +        ldiag/0/,idiag/40/,tdiag/ndiag*0./,
     +        dt/.025/,time/0./,gamma/1.667/,
     +        eta/0.02/,visx/0.02/,visy/0.02/,visz/0.02/
c uncomment the following lines for 
c   flux! computation or current! density output
c      data    tflux/0.0/,deltfl/1.0/,nflux/0/,
c     +        tjout/0.0/,deltj/1.0/,njout/0/
      end
	subroutine orgstart
c***********************************************************************
	include 'misflin'
c..........................................
c    startprozeduren:
c         -parameter read
c         -boundconditions
c         -initial state
c..........................................
      write(17,34)
      write(18,35)
      write(27,36)
      write(28,37)
      write(25,38)
c
      rewind 15
      read (15,1)  start,ieein,ieeout,jout,fluxcal,mitsat
      write(26,999)
  999 format(1x,'start read')
c    -----------------------------------------------------------
c     i) start of program by 
c            a: analytical state (start = true)
c            b: data from magtap (start = false)
c    -----------------------------------------------------------
      if(start)  then
	 call pread
	 iferror = 0
	 istep = 0
	 time = 0.
	 call grid
	 call bdcoef
	 call initcon
	 ferror = .false.
  150    continue
	 if ( (zentr(1)) .and. dt .lt. 0.5/difx((nx+1)/2) .or.
     +        (.not.(zentr(1))) .and. dt .lt. 0.5/difx(2) ) go to 200
	 dt = .5*dt
	 go to 150
  200    continue
      else
         call binin
	 call pread
	 call bdcoef
      end if
      if (nzz.gt.nz) nzz=nz
      call pwrite
      istep = 0
      
      startsat=.false.
      if (mitsat .and. tsat0 .eq. 0.) then 
       call satstart
       call satdat
      endif
c
      call bound
c uncomment the following lines for 
c   flux! computation 
c  or current! density output
c      if (fluxcal)  call initflux
c  grid for jout
c      do 910 ix = 1,nxj
c  910     xj(ix)=x(ix+ndiffj)
c      do 920 iy = 1,nyj
c  920     yj(iy)=y(iy+1)
c      do 930 iz = 1,nzj
c  930     zj(iz)=z(iz+1)
c      write(*,*) 'xj boundaries:', xj(1),xj(nxj)
c      write(*,*) 'yj boundaries:', yj(1),yj(nyj)
c      write(*,*) 'zj boundaries:', zj(1),zj(nzj)
  
c
    1 format(1x,l10)
    2 format(1x,i9)
    3 format(1x,f14.5)
c
   34 format(/1x,' time  !   masse    !  masse msp !',
     +           '  kin ener  !  mag ener  !  the energ',
     +       /1x,'-------!------------!------------!',
     +           '------------!------------!------------')
c
   35 format(/1x,' time  !v grad p/2!v dot jxb !res mal j2!',
     +                    'd/dt ekin ! d/dt eb  ! d/dt eth ',
     +       /1x,'-------!----------!----------!----------!',
     +                    '----------!----------!----------')
c
   36 format(/1x,' time  !   fluss    ! flussaend. !',
     +           ' aend ideal ! aend resist',
     +       /1x,'-------!------------!------------!',
     +           '------------!--------------')
c
   37 format(/1x,' time  !flux: mass! kin. en. ! the. en. !',
     +                    '  b. en.  !ideal.Con.! res. Con.',
     +       /1x,'-------!----------!----------!----------!',
     +                    '----------!----------!----------')
c
   38 format(/1x,' TIME  ',
     +       /1x,' MAXIMUM !   X   !    Y  !    Z  !',
     +           ' MINIMUM !   X   !    Y  !    Z   ',
     +       /1x,'---------!-------!-------!-------!',
     +           '---------!-------!-------!--------')
c
      return
      end
	subroutine intstart
c***********************************************************************
	include 'misflin'
c............................................................
c       integration der mhd-gleichungen:
c          1. integr. step = lax step on grid 1  
c          (igrid=1).  leap step for lax wendroff always on 
c          grid 0 for the leapfrog alternating on both grids
c............................................................
c
c    first two integrationssteps:
c...............................................
c
      write(26,61) istep
   61 format(1x,'output after initcon, integration step',i5)
c      call orgout
c      call symtx

      igrid = 1
      call lax
      call bound

      istep = istep + 1
      time = time + dt
      write(26,62) istep
   62 format(1x,'output after laxstep, integration step',i5)
c      call orgout
c      call symtx

      call termin
      if ( ferror .and. iferror.le.nsmooth ) call osmooth
      if ( ferror ) return
c
      igrid = 0
      call leap
      call bound
c
      istep = istep + 1
      time = time + dt
      write(26,63) istep
   63 format(1x,'output after leapstep, integration step',i5)
c      call orgout
c      call symtx

      call termin
      if ( ferror .and. iferror.le.nsmooth ) call osmooth
c
c      write(6,60) istep
c   60 format(1x,'integration step',i5)
c
      return
      end
	subroutine satstart
c***********************************************************************
	include 'misflin'
c
      integer    i,j,k,fnum,ier
      character*8      satname
c.................................................................
c    initializing of satparameters
c.................................................................
      startsat = .true.
      timesat = time
      open (31,file='satpar')
      
      do 50 i = 1,npos
       do 50 j = 1,nrow  
         k = j + (i-1)*nrow
         xsat(k) = xpos(i) + (j-1)*dxpos(i)
         ysat(k) = ypos(i) + (j-1)*dypos(i)
         zsat(k) = zpos(i) + (j-1)*dzpos(i)
         vxsat(k) = vxpos(i)
         vysat(k) = vypos(i)
         vzsat(k) = vzpos(i)
   50 continue

      write(31,21)  
      write(26,211)  
      write(31,22) nsat,timesat,tsatout
      write(26,22) nsat,timesat,tsatout
      write(31,23) nsim,bsim,vsim,psim,lsim,tsim
      write(31,24) rhobd(1),pbd(1),vxbd(1),vybd(1),vzbd(1),
     +             bxbd(1),bybd(1),bzbd(1),
     +             rhobd(2),pbd(2),vxbd(2),vybd(2),vzbd(2),
     +             bxbd(2),bybd(2),bzbd(2)
      write(31,25) 
      write(26,25) 
      do 60 k = 1,nsat
          write(31,26) k,xsat(k),ysat(k),zsat(k),
     +                     vxsat(k),vysat(k),vzsat(k)
          write(26,26) k,xsat(k),ysat(k),zsat(k),
     +                     vxsat(k),vysat(k),vzsat(k)
   60 continue
      write(31,27) 
      write(26,27) 
   
c  uncomment for cray
      if (ieein) then 
        call asnunit(49,'-a satbin -F f77 -N ieee',ier)
      else
        open (49,file='satbin',form='unformatted')
      end if
      write(49) nsat,timesat
      write(49) nsim,bsim,vsim,psim,lsim,tsim
      write(49) rhobd,pbd,vxbd,vybd,vzbd,bxbd,bybd,bzbd
      write(49) xsat,ysat,zsat,vxsat,vysat,vzsat

      do 70 i = 1,npos
      do 70 j = 1,nrow  
         k = j + (i-1)*nrow
         vxsat(k) = vxpos(i) + vxbd(1)
         vysat(k) = vypos(i) + vybd(1)
         vzsat(k) = vzpos(i) + vzbd(1)
   70 continue

   21 format(1x,'Satellite and program parameters:',
     +       /1x,'---------------------------------')
  211 format(1x,'Satellite parameters:',
     +       1x,'---------------------')
   22 format(/1x,'number of satellites:',i3,3x,
     +           'start time:',f7.2,3x,'time resolution:',f7.4)
   23 format(/1x,'Normalisation for : ',
     +       /1x,' No density ! magn. field!  velocity  !  pressure  !',
     +           'length units!    time     ',
     +       /1x,'  cm**(-3)  !     nT     !    km/s    !   nPascal  !',
     +           '     km     !      s      ',
     +         /,f11.2,f13.2,f13.2,f13.4,f13.2,f13.4)
   24 format(/1x,'Initial asymptotic values at xmin (1. row) ',
     +           'and xmax (2. row) for: ',
     +       /1x,'    rho      p        vx       vy       vz    ',
     +           '   bx       by       bz    ',
     +       /,8f9.3/,8f9.3)
   25 format(/1x,'Initial satellite locations and velocity:',
     +       /1x,'sat  !    x    !    y    !    z    !',
     +                 '    vx   !    vy   !    vz   ')
   26 format(1x,i3,'   ',f7.2,'   ',f7.2,'   ',f7.2,
     +             '   ',f7.2,'   ',f7.2,'   ',f7.2)
   27 format(1x,'Satellite velocity is relative to restframe of the',
     +            ' magnetosphere (at xmin)!!!',
     +       /1x,'Plasma velocity is recorded in the satellite frame!')

      return
      end
	subroutine satdat
c***********************************************************************
	include 'misflin'
c
	real    finterp
        integer ix,iy,iz,ixi,iyi,izi,i,j,k,ixsat,iysat,izsat,fnum
        real    dxsat(3),dysat(3),dzsat(3),
     +          bxint(2,2,2),byint(2,2,2),bzint(2,2,2), 
     +          sxint(2,2,2),syint(2,2,2),szint(2,2,2), 
     +          rhoint(2,2,2),pint(2,2,2),
     +          bxout,byout,bzout,sxout,syout,szout,rhoout,pout 
c.................................................................
      tsat0 = tsat0+tsatout
      do 5 k = 1,nsat
         xsat(k) = xsat(k) + (time-timesat)*vxsat(k)
         ysat(k) = ysat(k) + (time-timesat)*vysat(k)
         zsat(k) = zsat(k) + (time-timesat)*vzsat(k)
    5 continue

      do 50 k = 1,nsat
        ixsat = 2
        iysat = 2
        izsat = 2
        do 10 ix = 2,nx1
   10     if(xsat(k).ge.x(ix)) ixsat = ix
        do 20 iy = 2,ny1
   20     if(ysat(k).ge.y(iy)) iysat = iy
        do 30 iz = 2,nz1
   30     if(zsat(k).ge.z(iz)) izsat = iz
        dxsat(1) = x(ixsat+1)-xsat(k)
        dysat(1) = y(iysat+1)-ysat(k)
        dzsat(1) = z(izsat+1)-zsat(k)
        dxsat(2) = xsat(k)-x(ixsat)
        dysat(2) = ysat(k)-y(iysat)
        dzsat(2) = zsat(k)-z(izsat)
        dxsat(3) = x(ixsat+1)-x(ixsat)
        dysat(3) = y(iysat+1)-y(iysat)
        dzsat(3) = z(izsat+1)-z(izsat)

        if(xsat(k).lt.xmin) then
          xsat(k) = xmin
          ixsat = 2
          dxsat(1) = x(3)-x(2)
          dxsat(2) = 0.0
          dxsat(3) = x(3)-x(2)
        else if (xsat(k).gt.xmax) then
          xsat(k) = xmax
          ixsat = nx1
          dxsat(1) = x(nx)-x(nx1)
          dxsat(2) = 0.0
          dxsat(3) = x(nx)-x(nx1)
        end if

        if(ysat(k).lt.ymin) then
          if (perio(2)) then 
           ysat(k)=ysat(k)+(ymax-ymin)
          else
           ysat(k) = ymin
           iysat = 2
           dysat(1) = y(3)-y(2)
           dysat(2) = 0.0
           dysat(3) = y(3)-y(2)
          endif
        endif
        if (ysat(k).gt.ymax) then
          if (perio(2))  then
            ysat(k)=ysat(k)-(ymax-ymin)
          else
           ysat(k) = ymax
           iysat = ny1
           dysat(1) = y(ny)-y(ny1)
           dysat(2) = 0.0
           dysat(3) = y(ny)-y(ny1)
          endif
        end if

        if(zsat(k).lt.zmin) then
          if (perio(3)) then 
           zsat(k)=zsat(k)+(zmax-zmin)
          else
           zsat(k) = zmin
           izsat = 2
           dzsat(1) = z(3)-z(2)
           dzsat(2) = 0.0
           dzsat(3) = z(3)-z(2)
          endif
        endif
        if (zsat(k).gt.zmax) then
          if (perio(3))  then
            zsat(k)=zsat(k)-(zmax-zmin)
          else
          zsat(k) = zmax
          izsat = nz1
          dzsat(1) = z(nz)-z(nz1)
          dzsat(2) = 0.0
          dzsat(3) = z(nz)-z(nz1)
          endif
        end if 

        do 40 ixi = 0,1
        do 40 iyi = 0,1
        do 40 izi = 0,1
         bxint(ixi+1,iyi+1,izi+1)  = bx(ixi+ixsat,iyi+iysat,izi+izsat)
         byint(ixi+1,iyi+1,izi+1)  = by(ixi+ixsat,iyi+iysat,izi+izsat)
         bzint(ixi+1,iyi+1,izi+1)  = bz(ixi+ixsat,iyi+iysat,izi+izsat)
         sxint(ixi+1,iyi+1,izi+1)  = sx(ixi+ixsat,iyi+iysat,izi+izsat)
         syint(ixi+1,iyi+1,izi+1)  = sy(ixi+ixsat,iyi+iysat,izi+izsat)
         szint(ixi+1,iyi+1,izi+1)  = sz(ixi+ixsat,iyi+iysat,izi+izsat)
         rhoint(ixi+1,iyi+1,izi+1) = rho(ixi+ixsat,iyi+iysat,izi+izsat)
         pint(ixi+1,iyi+1,izi+1)   = 
     +                        2.*u(ixi+ixsat,iyi+iysat,izi+izsat)**gamma
   40   continue

        bxout = finterp(dxsat,dysat,dzsat,bxint)
        byout = finterp(dxsat,dysat,dzsat,byint)
        bzout = finterp(dxsat,dysat,dzsat,bzint)
        sxout = finterp(dxsat,dysat,dzsat,sxint)
        syout = finterp(dxsat,dysat,dzsat,syint)
        szout = finterp(dxsat,dysat,dzsat,szint)
        rhoout = finterp(dxsat,dysat,dzsat,rhoint)
        pout = finterp(dxsat,dysat,dzsat,pint)
	bxs(k) = bxout
	bys(k) = byout
	bzs(k) = bzout
	vxs(k) = sxout/rhoout
	vys(k) = syout/rhoout
	vzs(k) = szout/rhoout
	rhosat(k) = rhoout
	psat(k) = pout
	bsat(k) = sqrt( bxout*bxout+byout*byout+bzout*bzout )
	ptotsat(k) = psat(k) + bsat(k)*bsat(k)
   50 continue

      write(49) time,bxs,bys,bzs,vxs,vys,vzs,rhosat,psat,bsat,
     +               xsat,ysat,zsat
      
      timesat = time

   11 format(1x,f6.2,10f7.3)

      return
      end
      real function finterp(dx,dy,dz,fin)
c***********************************************************************
        dimension  dx(3),dy(3),dz(3),fin(2,2,2)
c.................................................................
      finterp = ( dx(1)*dy(1)*dz(1)*fin(1,1,1)
     +       + dx(1)*dy(1)*dz(2)*fin(1,1,2)
     +       + dx(1)*dy(2)*dz(1)*fin(1,2,1)
     +       + dx(1)*dy(2)*dz(2)*fin(1,2,2)
     +       + dx(2)*dy(1)*dz(1)*fin(2,1,1)
     +       + dx(2)*dy(1)*dz(2)*fin(2,1,2)
     +       + dx(2)*dy(2)*dz(1)*fin(2,2,1)
     +       + dx(2)*dy(2)*dz(2)*fin(2,2,2) )
     +                          / ( dx(3)*dy(3)*dz(3) )
      return
      end
	subroutine orgint
c***********************************************************************
	include 'misflin'
c................................................................
c
c  1)  leapfrog integration:
c.................................
c
      if (intart .eq. 1) then
	 igrid = 1
	 call leap
	 call bound
	 istep = istep + 1
	 time = time + dt
c	 call symtx
	 call termin
	 if ( ferror .and. iferror.le.nsmooth ) call osmooth
	 if ( ferror ) return

	 igrid = 0
	 call leap
	 call bound
	 istep = istep + 1
	 time = time + dt
	 if (mod(istep,100).eq.0) write(6,60) istep
c	 call symtx
	 call termin
	 if ( ferror .and. iferror.le.nsmooth ) call osmooth
      end if
c
c
c  2)  lax-wendroff integration:
c.....................................
c
      if (intart .eq. 2) then
	 igrid = 1
	 call meanlax
	 call bound
c	 call symtx
	 call lax
	 call bound
	 istep = istep + 1
	 time = time + dt
c	 call symtx
	 call termin
	 if ( ferror .and. iferror.le.nsmooth ) call osmooth
	 if ( ferror ) return

	 igrid = 0
	 call leap
	 call bound
	 istep = istep + 1
	 time = time + dt
c	 write(6,60) istep
c 	 call symtx
	 call termin
	 if ( ferror .and. iferror.le.nsmooth ) call osmooth
      end if
c
   60 format(1x,'integration step',i5)
c
      return
      end
	subroutine osmooth
c***********************************************************************
	include 'misflin'
c
c.................................................................
      iferror = iferror + 1
      write(26,1) iferror
c      write(6,1) iferror
      ferror = .false.
      call smooth
      call bound
      call termin
c
    1 format(1x,'anzahl der glaettungen: ',i3)
      return
      end
	subroutine pread
c***********************************************************************
	include 'misflin'
c
	logical  ldum
	integer  idum, i, j
	real     rdum
c.....................................................
      read (15,2) etasw,isafe,iend,intart,ivisrho,ivisu,nsmooth,idiag
      write (6,2) etasw,isafe,iend,intart,ivisrho,ivisu,nsmooth,idiag
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
	subroutine pwrite
c***********************************************************************
	include 'misflin'
c
	integer  i, j, ix, iy, iz
c.....................................................
      write(26,10)  istep
      write(26,11)  nx,ny,nz,
     +               intart,isafe,iend,nsmooth,ivisrho,ivisu,
     +               igrid,idiag,iferror,
     +               start,ieein,ieeout,jout,fluxcal,mitsat,deltj
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
	subroutine grid
c***********************************************************************
	include 'misflin'
c
	integer  ix,ixs,iy,iys,iz,ix0,iy0,iz0
	real     kx,ky,kz,xanf,yanf,zanf,
     +           beta,a1,a2,a3,a4,w,w1,w2
c .................................................................
c     grid computes gridparameters
c         ( x perp to current sheet       -  n koordinate
c           y                             - -m koordinate
c           z anti parallel field         -  l koordinate )
c .................................................................
c
c   x-coordinate
c  (y und z similar) using:
c         x  =  w + a1*sin(beta*w) + a2*sin(2*beta*w)
c         w  =  kx*(ix-ix0)
c         kx =  2.*xmax/(nx-3) ; ix0 = (nx+1)/2
c         beta = pi/xmax
c         a1,a2 determined by grid distance at x=0 (eps(1))
c         difx - coeff 1. derivative
c         ddifx, ddifpx und ddifmx - 2. derivative
c ..............................................................
      write(*,991)
 991  format('vor x gitter')
      kx =  (xmax-xmin)/(nx-3)
      xanf = xmin
      ix0 = 2
      beta = pi/(xmax-xmin)
      if (zentr(1)) then
	xanf = (xmin+xmax)/2.0
	ix0 = (nx+1)/2
	beta = 2*pi/(xmax-xmin)
      end if
      a1 = 0.0
      a2 = 0.0
      a3 = 0.0
      a4 = 0.0
      if (.not. aequi(1)) then
	a1 = 24./17.*(eps(1)-kx)/kx/beta
	a2 =  6./17.*(kx-eps(1))/kx/beta
	a3 =  8./51.*(eps(1)-kx)/kx/beta
	a4 =  3./68.*(kx-eps(1))/kx/beta
      end if
c.......................................
      do 200 ix = 1,nx
         w      =  kx*(ix-ix0)
	 x(ix)  =  w + a1*sin(beta*w) + a2*sin(2*beta*w)
     +               + a3*sin(3.*beta*w) + a4*sin(4.*beta*w)
	 difx(ix) = 0.5/kx/(1. + a1*beta*cos(beta*w)
     +                    + 2.*a2*beta*cos(2.*beta*w)
     +                     + 3.*a3*beta*cos(3.*beta*w)
     +                      + 4.*a4*beta*cos(4.*beta*w) )
         w1 = w + 0.5*kx
	 w2 = w - 0.5*kx
	 ddifpx(ix) = 2.*difx(ix)/kx/(1. + a1*beta*cos(beta*w1)
     +                              + 2.*a2*beta*cos(2.*beta*w1)
     +                               + 3.*a3*beta*cos(3.*beta*w1)
     +                                + 4.*a4*beta*cos(4.*beta*w1) )
	 ddifmx(ix) = 2.*difx(ix)/kx/(1. + a1*beta*cos(beta*w2)
     +                              + 2.*a2*beta*cos(2.*beta*w2)
     +                               + 3.*a3*beta*cos(3.*beta*w2)
     +                                + 4.*a4*beta*cos(4.*beta*w2) )
         ddifx(ix) = .5*(ddifpx(ix)+ddifmx(ix))
  200 continue
c
c  coeff for averages:
c...........................................
      do 250 ix = 2,nx1
	 meanmx(ix) = (x(ix+1)-x(ix))/(x(ix+1)-x(ix-1))
	 meanpx(ix) = (x(ix)-x(ix-1))/(x(ix+1)-x(ix-1))
  250 continue
      meanmx(1)   = meanpx(3)
      meanpx(1)   = meanmx(3)
      meanmx(nx)  = meanpx(nx-2)
      meanpx(nx)  = meanmx(nx-2)
c
c symmetriy for centered nonuniform grid
c.......................................................
      if (zentr(1))   then
	do 260 ix = 1, ix0
	  x(ix)      = ( x(ix)-x(nx+1-ix) )/2.0
	  difx(ix)   = ( difx(ix)+difx(nx+1-ix) )/2.0
	  ddifpx(ix) = ( ddifpx(ix)+ddifmx(nx+1-ix) )/2.0
	  ddifmx(ix) = ( ddifmx(ix)+ddifpx(nx+1-ix) )/2.0
	  ddifx(ix)  = ( ddifx(ix)+ddifx(nx+1-ix) )/2.0
	  meanpx(ix) = ( meanpx(ix)+meanmx(nx+1-ix) )/2.0
	  meanmx(ix) = ( meanmx(ix)+meanpx(nx+1-ix) )/2.0
  260   continue
	  x(ix0) = 0.0
	  ddifpx(ix0) = ( ddifpx(ix0)+ddifmx(ix0) )/2.0
	  ddifmx(ix0) = ( ddifpx(ix0)+ddifmx(ix0) )/2.0
	  ddifx(ix0) =  ( ddifpx(ix0)+ddifmx(ix0) )/2.0
	  meanpx(ix0) = ( meanpx(ix0)+meanmx(ix0) )/2.0
	  meanmx(ix0) = ( meanpx(ix0)+meanmx(ix0) )/2.0
	do 270 ix = ix0+1, nx
	  x(ix)      = -x(nx+1-ix)
	  difx(ix)   = difx(nx+1-ix)
	  ddifpx(ix) = ddifmx(nx+1-ix)
	  ddifmx(ix) = ddifpx(nx+1-ix)
	  ddifx(ix)  = ddifx(nx+1-ix)
	  meanpx(ix) = meanmx(nx+1-ix)
	  meanmx(ix) = meanpx(nx+1-ix)
  270   continue

      end if
      do 300 ix = 1, nx
	x(ix) = xanf + x(ix)
  300 continue
c
c  gridboundary
c..............................................
      x(1)        = 2.*xmin - x(3)
      x(nx)       = 2.*xmax -x(nx-2)
      x(2)        = xmin
      x(nx1)     = xmax
      difx(1)     = difx(3)
      difx(nx)    = difx(nx-2)
      ddifmx(2)   = ( ddifpx(2) + ddifmx(2) ) / 2.0
      ddifmx(nx1)= ( ddifpx(nx1) + ddifmx(nx1) ) / 2.0
      ddifpx(2)   = ddifmx(2)
      ddifpx(nx1)= ddifmx(nx1)
      ddifmx(1)   = ddifpx(3)
      ddifpx(1)   = ddifmx(3)
      ddifmx(nx)  = ddifpx(nx-2)
      ddifpx(nx)  = ddifmx(nx-2)
      meanmx(1)   = meanpx(3)
      meanpx(1)   = meanmx(3)
      meanmx(nx)  = meanpx(nx-2)
      meanpx(nx)  = meanmx(nx-2)
c
      if (zentr(1))   then
	do 360 ix = 1,nx
	  ixs = 1+nx-ix
	  if ( meanpx(ix) .ne. meanmx(ixs) )
     +      write(26,21) ix,meanpx(ix),meanmx(ixs),time
	  if ( ddifpx(ix) .ne. ddifmx(ixs) )
     +      write(26,22) ix,ddifpx(ix),ddifmx(ixs),time
  360   continue
      end if
   21 format(' meanpx(',i3,') =',f12.7,
     +       '   and',f12.7,'   time =',f9.4)
   22 format(' ddifpx(',i3,') =',f12.7,
     +       '   and',f12.7,'   time =',f9.4)

c.......................................................
c   y-coordinate
c      analogue x
c.......................................................
      ky =  (ymax-ymin)/(ny-3)
      yanf = ymin
      iy0 = 2
      beta = pi/(ymax-ymin)
      if (zentr(2)) then
	yanf = (ymin+ymax)/2.0
	iy0 = (ny+1)/2
	beta = 2*pi/(ymax-ymin)
      end if
      a1 = 0.0
      a2 = 0.0
      a3 = 0.0
      a4 = 0.0
      if (.not. aequi(2)) then
	a1 = 24./17.*(eps(2)-ky)/ky/beta
	a2 =  6./17.*(ky-eps(2))/ky/beta
	a3 =  8./51.*(eps(2)-ky)/ky/beta
	a4 =  3./68.*(ky-eps(2))/ky/beta
      end if
c.......................................
      do 400 iy = 1,ny
         w      =  ky*(iy-iy0)
	 y(iy)  =  w + a1*sin(beta*w) + a2*sin(2*beta*w)
     +               + a3*sin(3.*beta*w) + a4*sin(4.*beta*w)
	 dify(iy) = 0.5/ky/(1. + a1*beta*cos(beta*w)
     +                    + 2.*a2*beta*cos(2.*beta*w)
     +                     + 3.*a3*beta*cos(3.*beta*w)
     +                      + 4.*a4*beta*cos(4.*beta*w) )
         w1 = w + 0.5*ky
         w2 = w - 0.5*ky
	 ddifpy(iy) = 2.*dify(iy)/ky/(1. + a1*beta*cos(beta*w1)
     +                              + 2.*a2*beta*cos(2.*beta*w1)
     +                               + 3.*a3*beta*cos(3.*beta*w1)
     +                                + 4.*a4*beta*cos(4.*beta*w1) )
	 ddifmy(iy) = 2.*dify(iy)/ky/(1. + a1*beta*cos(beta*w2)
     +                              + 2.*a2*beta*cos(2.*beta*w2)
     +                               + 3.*a3*beta*cos(3.*beta*w2)
     +                                + 4.*a4*beta*cos(4.*beta*w2) )
	 ddify(iy) = .5*(ddifpy(iy)+ddifmy(iy))
  400 continue
c
c  coeff for averages
c............................................
c
      do 450 iy = 2,ny1
	 meanmy(iy) = (y(iy+1)-y(iy))/(y(iy+1)-y(iy-1))
	 meanpy(iy) = (y(iy)-y(iy-1))/(y(iy+1)-y(iy-1))
  450 continue
      meanmy(1) = meanpy(3)
      meanpy(1) = meanmy(3)
      meanmy(ny)   = meanpy(ny-2)
      meanpy(ny)   = meanmy(ny-2)
c
c centered nonuniform grid
c.......................................................
      if (zentr(2))   then
	do 460 iy = 1, iy0
	  y(iy)      = ( y(iy)-y(ny+1-iy) )/2.0
	  dify(iy)   = ( dify(iy)+dify(ny+1-iy) )/2.0
	  ddifpy(iy) = ( ddifpy(iy)+ddifmy(ny+1-iy) )/2.0
	  ddifmy(iy) = ( ddifmy(iy)+ddifpy(ny+1-iy) )/2.0
	  ddify(iy)  = ( ddify(iy)+ddify(ny+1-iy) )/2.0
	  meanpy(iy) = ( meanpy(iy)+meanmy(ny+1-iy) )/2.0
	  meanmy(iy) = ( meanmy(iy)+meanpy(ny+1-iy) )/2.0
  460   continue
	  y(iy0) = 0.0
	  ddifpy(iy0) = ( ddifpy(iy0)+ddifmy(iy0) )/2.0
	  ddifmy(iy0) = ( ddifpy(iy0)+ddifmy(iy0) )/2.0
	  ddify(iy0) =  ( ddifpy(iy0)+ddifmy(iy0) )/2.0
	  meanpy(iy0) = ( meanpy(iy0)+meanmy(iy0) )/2.0
	  meanmy(iy0) = ( meanpy(iy0)+meanmy(iy0) )/2.0
	do 470 iy = iy0+1, ny
	  y(iy)      = -y(ny+1-iy)
	  dify(iy)   = dify(ny+1-iy)
	  ddifpy(iy) = ddifmy(ny+1-iy)
	  ddifmy(iy) = ddifpy(ny+1-iy)
	  ddify(iy)  = ddify(ny+1-iy)
	  meanpy(iy) = meanmy(ny+1-iy)
	  meanmy(iy) = meanpy(ny+1-iy)
  470   continue
      end if
      do 500 iy = 1, ny
	y(iy) = yanf + y(iy)
  500 continue
c
c  grid boundary
c...............................................
      y(1)        = 2.*ymin - y(3)
      y(ny)       = 2.*ymax - y(ny-2)
      y(2)        = ymin
      y(ny1)     = ymax
      dify(1)     = dify(3)
      dify(ny)    = dify(ny-2)
      ddifmy(2)   = ( ddifpy(2) + ddifmy(2) ) / 2.0
      ddifmy(ny1)= ( ddifpy(ny1) + ddifmy(ny1) ) / 2.0
      ddifpy(2)   = ddifmy(2)
      ddifpy(ny1)= ddifmy(ny1)
      ddifmy(1)   = ddifpy(3)
      ddifpy(1)   = ddifmy(3)
      ddifmy(ny)  = ddifpy(ny-2)
      ddifpy(ny)  = ddifmy(ny-2)
      meanmy(1)   = meanpy(3)
      meanpy(1)   = meanmy(3)
      meanmy(ny)  = meanpy(ny-2)
      meanpy(ny)  = meanmy(ny-2)
c
      if (zentr(2))   then
	do 560 iy = 1,ny
	  iys = 1+ny-iy
	  if ( meanpy(iy) .ne. meanmy(iys) )
     +      write(26,31) iy,meanpy(iy),meanmy(iys),time
	  if ( ddifpy(iy) .ne. ddifmy(iys) )
     +      write(26,32) iy,ddifpy(iy),ddifmy(iys),time
  560   continue
      end if
   31 format(' meanpy(',i3,') =',f12.7,
     +       '   und',f12.7,'   time =',f9.4)
   32 format(' ddifpy(',i3,') =',f12.7,
     +       '   und',f12.7,'   time =',f9.4)

c.......................................................
c  z-coordinate
c      analogue  x
c.......................................................
      kz =  (zmax-zmin)/(nz-3)
      zanf = zmin
      iz0 = 2
      beta = pi/(zmax-zmin)
      if (zentr(3)) then
	zanf = (zmin+zmax)/2.0
	iz0 = (nz+1)/2
	beta = 2*pi/(zmax-zmin)
      end if
      a1 = 0.0
      a2 = 0.0
      a3 = 0.0
      a4 = 0.0
      if (.not. aequi(3)) then
	a1 = 24./17.*(eps(3)-kz)/kz/beta
	a2 =  6./17.*(kz-eps(3))/kz/beta
	a3 =  8./51.*(eps(3)-kz)/kz/beta
	a4 =  3./68.*(kz-eps(3))/kz/beta
      end if
c.......................................
      do 600 iz = 1,nz
         w      =  kz*(iz-iz0)
	 z(iz)  =  w + a1*sin(beta*w) + a2*sin(2*beta*w)
     +               + a3*sin(3.*beta*w) + a4*sin(4.*beta*w)
	 difz(iz) = 0.5/kz/(1. + a1*beta*cos(beta*w)
     +                    + 2.*a2*beta*cos(2.*beta*w)
     +                     + 3.*a3*beta*cos(3.*beta*w)
     +                      + 4.*a4*beta*cos(4.*beta*w) )
         w1 = w + 0.5*kz
         w2 = w - 0.5*kz
	 ddifpz(iz) = 2.*difz(iz)/kz/(1. + a1*beta*cos(beta*w1)
     +                              + 2.*a2*beta*cos(2.*beta*w1)
     +                               + 3.*a3*beta*cos(3.*beta*w1)
     +                                + 4.*a4*beta*cos(4.*beta*w1) )
	 ddifmz(iz) = 2.*difz(iz)/kz/(1. + a1*beta*cos(beta*w2)
     +                              + 2.*a2*beta*cos(2.*beta*w2)
     +                               + 3.*a3*beta*cos(3.*beta*w2)
     +                                + 4.*a4*beta*cos(4.*beta*w2) )
         ddifz(iz) = .5*(ddifpz(iz)+ddifmz(iz))
  600 continue
c
c  averages
c............................................
      do 650 iz = 2,nz1
	 meanmz(iz) = (z(iz+1)-z(iz))/(z(iz+1)-z(iz-1))
         meanpz(iz) = (z(iz)-z(iz-1))/(z(iz+1)-z(iz-1))
  650 continue
      meanmz(1)   = meanpz(3)
      meanpz(1)   = meanmz(3)
      meanmz(nz)  = meanpz(nz-2)
      meanpz(nz)  = meanmz(nz-2)
c
c centered nonuniform grid
c.......................................................
      if (zentr(3))   then
	do 660 iz = 1, iz0
	  z(iz)      = ( z(iz)-z(nz+1-iz) )/2.0
	  difz(iz)   = ( difz(iz)+difz(nz+1-iz) )/2.0
	  ddifpz(iz) = ( ddifpz(iz)+ddifmz(nz+1-iz) )/2.0
	  ddifmz(iz) = ( ddifmz(iz)+ddifpz(nz+1-iz) )/2.0
	  ddifz(iz)  = ( ddifz(iz)+ddifz(nz+1-iz) )/2.0
	  meanpz(iz) = ( meanpz(iz)+meanmz(nz+1-iz) )/2.0
	  meanmz(iz) = ( meanmz(iz)+meanpz(nz+1-iz) )/2.0
  660   continue
	  z(iz0) = 0.0
	  ddifpz(iz0) = ( ddifpz(iz0)+ddifmz(iz0) )/2.0
	  ddifmz(iz0) = ( ddifpz(iz0)+ddifmz(iz0) )/2.0
	  ddifz(iz0) =  ( ddifpz(iz0)+ddifmz(iz0) )/2.0
	  meanpz(iz0) = ( meanpz(iz0)+meanmz(iz0) )/2.0
	  meanmz(iz0) = ( meanpz(iz0)+meanmz(iz0) )/2.0
	do 670 iz = iz0+1, nz
	  z(iz)      = -z(nz+1-iz)
	  difz(iz)   = difz(nz+1-iz)
	  ddifpz(iz) = ddifmz(nz+1-iz)
	  ddifmz(iz) = ddifpz(nz+1-iz)
	  ddifz(iz)  = ddifz(nz+1-iz)
	  meanpz(iz) = meanmz(nz+1-iz)
	  meanmz(iz) = meanpz(nz+1-iz)
  670   continue
      end if
      do 700 iz = 1, nz
	z(iz) = zanf + z(iz)
  700 continue
c
c  grid boundary
c..............................................
      z(1)        = 2.*zmin - z(3)
      z(nz)       = 2.*zmax -z(nz-2)
      z(2)        = zmin
      z(nz1)     = zmax
      difz(1)     = difz(3)
      difz(nz)    = difz(nz-2)
      ddifmz(2)   = ( ddifpz(2) + ddifmz(2) ) / 2.0
      ddifmz(nz1)= ( ddifpz(nz1) + ddifmz(nz1) ) / 2.0
      ddifpz(2)   = ddifmz(2)
      ddifpz(nz1)= ddifmz(nz1)
      ddifmz(1)   = ddifpz(3)
      ddifpz(1)   = ddifmz(3)
      ddifmz(nz)  = ddifpz(nz-2)
      ddifpz(nz)  = ddifmz(nz-2)
      meanmz(1)   = meanpz(3)
      meanpz(1)   = meanmz(3)
      meanmz(nz)  = meanpz(nz-2)
      meanpz(nz)  = meanmz(nz-2)

      if (zentr(3)) then
	do 760 iz = 1,nz
	  if ( meanpz(iz) .ne. meanmz(nz+1-iz) )
     +      write(26,41) iz,meanpz(iz),meanmz(nz+1-iz),time
	  if ( ddifpz(iz) .ne. ddifmz(nz+1-iz) )
     +      write(26,42) iz,ddifpz(iz),ddifmz(nz+1-iz),time
  760   continue
      end if
   41 format(' meanpz(',i3,') =',f12.7,
     +       '   und',f12.7,'   time =',f9.4)
   42 format(' ddifpz(',i3,') =',f12.7,
     +       '   und',f12.7,'   time =',f9.4)
c
      return
      end
	subroutine initcon
c***********************************************************************
	include 'misflin'
c
	integer  ix,iy,iz
	real     gaminv,phirad,psirad,eta0,sinpsi,cospsi,
     +           pprof(nx)
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
	 syprof(ix)  = 0.0
	 syprof(ix)  = 0.0
   50 continue
      do 55 ix = 1,nx
	 uprof(ix) = (pprof(ix)/2.)**gaminv
   55 continue
      do 60 iz = 1,nz
	do 60 iy = 1,ny
	do 60 ix = 1,nx
	  rho(ix,iy,iz) = rhoprof(ix)
	  help(ix,iy,iz)= pprof(ix)
	  bx(ix,iy,iz)  = bxprof(ix)
	  by(ix,iy,iz)  = byprof(ix)
	  bz(ix,iy,iz)  = bzprof(ix)
	  sx(ix,iy,iz)  = sxprof(ix)
	  sy(ix,iy,iz)  = syprof(ix)
	  sz(ix,iy,iz)  = szprof(ix)
	  prof(ix,iy,iz) = (1. - tanh(x(ix))**2)
     +                    *(1. - tanh((z(iz)-0.)/3.)**2)
c     +       *0.5*( tanh((y(iy)+10.)/4.) - tanh((y(iy)-10.)/4.) )
     +                    *(1. - tanh(y(iy)/5.)**2)
   60 continue
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
	    hilfy(ix,iy,iz) = sy(ix,iy,iz)
	    hilfz(ix,iy,iz) = sz(ix,iy,iz)
	    feldy(ix,iy,iz) = by(ix,iy,iz)
	    feldz(ix,iy,iz) = bz(ix,iy,iz)
  110  continue
       do 120 iz = 1,nz
       do 120 iy = 1,ny
        do 120 ix = 1,nx
	 sy(ix,iy,iz) = cospsi*hilfy(ix,iy,iz) + sinpsi*hilfz(ix,iy,iz)
	 sz(ix,iy,iz) = -sinpsi*hilfy(ix,iy,iz) + cospsi*hilfz(ix,iy,iz)
	 by(ix,iy,iz) = cospsi*feldy(ix,iy,iz) + sinpsi*feldz(ix,iy,iz)
	 bz(ix,iy,iz) = -sinpsi*feldy(ix,iy,iz) + cospsi*feldz(ix,iy,iz)
  120  continue
      endif
c
c    boundary values at xmin and xmax
c
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
	subroutine smooth
c***********************************************************************
	include 'misflin'
c
	integer  ix,iy,iz
	real     ppp1(nx),ppp2(nx)
c ................................................................
	do 50 ix = 1,nx
	  ppp1(ix) =  (1.-0.75/cosh(4.*x(ix)/xmax))
          ppp2(ix) =  0.25/cosh(4.*x(ix)/xmax)
   50   continue
      do 100 iz = 2,nz1
	do 100 iy = 2,ny1
	do 100 ix = 2,nx1
	  help(ix,iy,iz) =  ppp1(ix)*rho(ix,iy,iz)
     +              + ppp2(ix)*(meanpx(ix)*rho(ix+1,iy,iz)
     +                        + meanmx(ix)*rho(ix-1,iy,iz)
     +                        + meanpy(iy)*rho(ix,iy+1,iz)
     +                        + meanmy(iy)*rho(ix,iy-1,iz)
     +                        + meanpz(iz)*rho(ix,iy,iz+1)
     +                        + meanmz(iz)*rho(ix,iy,iz-1) )
  100 continue
      do 101 iz = 2,nz1
	do 101 iy = 2,ny1
	do 101 ix = 2,nx1
	    rho(ix,iy,iz) = help(ix,iy,iz)
  101 continue
c
      do 200 iz = 2,nz1
	do 200 iy = 2,ny1
	do 200 ix = 2,nx1
	  help(ix,iy,iz) =  ppp1(ix)*u(ix,iy,iz)
     +              + ppp2(ix)*(meanpx(ix)*u(ix+1,iy,iz)
     +                        + meanmx(ix)*u(ix-1,iy,iz)
     +                        + meanpy(iy)*u(ix,iy+1,iz)
     +                        + meanmy(iy)*u(ix,iy-1,iz)
     +                        + meanpz(iz)*u(ix,iy,iz+1)
     +                        + meanmz(iz)*u(ix,iy,iz-1) )
  200 continue
      do 201 iz = 2,nz1
	do 201 iy = 2,ny1
	do 201 ix = 2,nx1
	    u(ix,iy,iz) = help(ix,iy,iz)
  201 continue
c
      do 300 iz = 2,nz1
	do 300 iy = 2,ny1
	do 300 ix = 2,nx1
	  feldx(ix,iy,iz) =  ppp1(ix)*sx(ix,iy,iz)
     +              + ppp2(ix)*(meanpx(ix)*sx(ix+1,iy,iz)
     +                        + meanmx(ix)*sx(ix-1,iy,iz)
     +                        + meanpy(iy)*sx(ix,iy+1,iz)
     +                        + meanmy(iy)*sx(ix,iy-1,iz)
     +                        + meanpz(iz)*sx(ix,iy,iz+1)
     +                        + meanmz(iz)*sx(ix,iy,iz-1) )
	  feldy(ix,iy,iz) =  ppp1(ix)*sy(ix,iy,iz)
     +              + ppp2(ix)*(meanpx(ix)*sy(ix+1,iy,iz)
     +                        + meanmx(ix)*sy(ix-1,iy,iz)
     +                        + meanpy(iy)*sy(ix,iy+1,iz)
     +                        + meanmy(iy)*sy(ix,iy-1,iz)
     +                        + meanpz(iz)*sy(ix,iy,iz+1)
     +                        + meanmz(iz)*sy(ix,iy,iz-1) )
	  feldz(ix,iy,iz) =  ppp1(ix)*sz(ix,iy,iz)
     +              + ppp2(ix)*(meanpx(ix)*sz(ix+1,iy,iz)
     +                        + meanmx(ix)*sz(ix-1,iy,iz)
     +                        + meanpy(iy)*sz(ix,iy+1,iz)
     +                        + meanmy(iy)*sz(ix,iy-1,iz)
     +                        + meanpz(iz)*sz(ix,iy,iz+1)
     +                        + meanmz(iz)*sz(ix,iy,iz-1) )
  300 continue
      do 301 iz = 2,nz1
	do 301 iy = 2,ny1
	do 301 ix = 2,nx1
	    sx(ix,iy,iz) = feldx(ix,iy,iz)
	    sy(ix,iy,iz) = feldy(ix,iy,iz)
	    sz(ix,iy,iz) = feldz(ix,iy,iz)
  301 continue
c
      do 400 iz = 2,nz1
	do 400 iy = 2,ny1
	do 400 ix = 2,nx1
	  feldx(ix,iy,iz) =  ppp1(ix)*bx(ix,iy,iz)
     +              + ppp2(ix)*(meanpx(ix)*bx(ix+1,iy,iz)
     +                        + meanmx(ix)*bx(ix-1,iy,iz)
     +                        + meanpy(iy)*bx(ix,iy+1,iz)
     +                        + meanmy(iy)*bx(ix,iy-1,iz)
     +                        + meanpz(iz)*bx(ix,iy,iz+1)
     +                        + meanmz(iz)*bx(ix,iy,iz-1) )
	  feldy(ix,iy,iz) =  ppp1(ix)*by(ix,iy,iz)
     +              + ppp2(ix)*(meanpx(ix)*by(ix+1,iy,iz)
     +                        + meanmx(ix)*by(ix-1,iy,iz)
     +                        + meanpy(iy)*by(ix,iy+1,iz)
     +                        + meanmy(iy)*by(ix,iy-1,iz)
     +                        + meanpz(iz)*by(ix,iy,iz+1)
     +                        + meanmz(iz)*by(ix,iy,iz-1) )
	  feldz(ix,iy,iz) =  ppp1(ix)*bz(ix,iy,iz)
     +              + ppp2(ix)*(meanpx(ix)*bz(ix+1,iy,iz)
     +                        + meanmx(ix)*bz(ix-1,iy,iz)
     +                        + meanpy(iy)*bz(ix,iy+1,iz)
     +                        + meanmy(iy)*bz(ix,iy-1,iz)
     +                        + meanpz(iz)*bz(ix,iy,iz+1)
     +                        + meanmz(iz)*bz(ix,iy,iz-1) )
  400 continue
      do 401 iz = 2,nz1
	do 401 iy = 2,ny1
	do 401 ix = 2,nx1
	    bx(ix,iy,iz) = feldx(ix,iy,iz)
	    by(ix,iy,iz) = feldy(ix,iy,iz)
	    bz(ix,iy,iz) = feldz(ix,iy,iz)
  401 continue
c
c
      return
      end
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
	subroutine lax
c***********************************************************************
	include 'misflin'
c
	integer  ix,iy,iz,ixanf
	real     eingam,chelp
c .........................................................
c  integration of 3-d mhd equations
c .........................................................
      eingam = 1 - gamma
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
     +            + difz(iz)*(feldy(ix,iy,iz+1)-feldy(ix,iy,iz-1)) )
	  sy(ix,iy,iz) = sy(ix,iy,iz) - dt*(
     +              dify(iy)*(hilfy(ix,iy+1,iz)-hilfy(ix,iy-1,iz))
     +            + difz(iz)*(feldx(ix,iy,iz+1)-feldx(ix,iy,iz-1))
     +            + difx(ix)*(feldz(ix+1,iy,iz)-feldz(ix-1,iy,iz)) )
	  sz(ix,iy,iz) = sz(ix,iy,iz) - dt*(
     +              difz(iz)*(hilfz(ix,iy,iz+1)-hilfz(ix,iy,iz-1))
     +            + difx(ix)*(feldy(ix+1,iy,iz)-feldy(ix-1,iy,iz))
     +            + dify(iy)*(feldx(ix,iy+1,iz)-feldx(ix,iy-1,iz)) )
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
c
      return
      end
	subroutine leap
c***********************************************************************
	include 'misflin'
c
	integer  ix,iy,iz,ixanf
	real     ddt,eingam,dif0,chelp
c .........................................................
c  integration of 3-d mhdequations
c .........................................................
      ddt = 2.*dt
      eingam = 1 - gamma
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
      do 245 iz = 2,nzz-1
      do 245 iy = 2,ny1
      ixanf = 2 + mod(iz+iy+igrid,2)
      do 245 ix = ixanf,nx1,2
	 hilfx(ix,iy,iz)=difx(ix)*(hilfx(ix+1,iy,iz)-hilfx(ix-1,iy,iz))
     +                 + dify(iy)*(feldz(ix,iy+1,iz)-feldz(ix,iy-1,iz))
     +                 + difz(iz)*(feldy(ix,iy,iz+1)-feldy(ix,iy,iz-1))
	 hilfy(ix,iy,iz)=dify(iy)*(hilfy(ix,iy+1,iz)-hilfy(ix,iy-1,iz))
     +                 + difz(iz)*(feldx(ix,iy,iz+1)-feldx(ix,iy,iz-1))
     +                 + difx(ix)*(feldz(ix+1,iy,iz)-feldz(ix-1,iy,iz))
         hilfz(ix,iy,iz)=difz(iz)*(hilfz(ix,iy,iz+1)-hilfz(ix,iy,iz-1))
     +                 + difx(ix)*(feldy(ix+1,iy,iz)-feldy(ix-1,iy,iz))
     +                 + dify(iy)*(feldx(ix,iy+1,iz)-feldx(ix,iy-1,iz))
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
c      hilf = v x b
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
      do 830 iz = 2,nzz-1
      do 830 iy = 2,ny1
      ixanf = 2 + mod(iz+iy+igrid,2)
      do 830 ix = ixanf,nx1,2
	 help(ix,iy,iz) = difx(ix)*(hilfx(ix+1,iy,iz)-hilfx(ix-1,iy,iz))
     +         + dify(iy)*(hilfy(ix,iy+1,iz)-hilfy(ix,iy-1,iz))
     +         + difz(iz)*(hilfz(ix,iy,iz+1)-hilfz(ix,iy,iz-1))
     +         + chelp*res(ix,iy,iz)*hilf(ix,iy,iz)*(
     +             ( meanpx(ix)*u(ix+1,iy,iz)+meanmx(ix)*u(ix-1,iy,iz) )
     +            +( meanpy(iy)*u(ix,iy+1,iz)+meanmy(iy)*u(ix,iy-1,iz) )
     +            +( meanpz(iz)*u(ix,iy,iz+1)+meanmz(iz)*u(ix,iy,iz-1) )
     +                           )**eingam 

  830 continue
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
c
      return
      end
	subroutine fbound0
c***********************************************************************
	include 'misflin'
c
	integer  ix,iy,iz
c .........................................................
c  set boundary of feld to 0
c .........................................................
	 do 121 iy = 1,ny
	 do 121 ix = 1,nx
	     feldx(ix,iy,1) = 0.
	     feldy(ix,iy,1) = 0.
	     feldz(ix,iy,1) = 0.
	     feldx(ix,iy,nz) = 0.
	     feldy(ix,iy,nz) = 0.
	     feldz(ix,iy,nz) = 0.
  121    continue
	 do 122 iz = 1,nz
	   do 122 ix = 1,nx
	     feldx(ix,1,iz) = 0.
	     feldy(ix,1,iz) = 0.
	     feldz(ix,1,iz) = 0.
	     feldx(ix,ny,iz) = 0.
	     feldy(ix,ny,iz) = 0.
	     feldz(ix,ny,iz) = 0.
  122    continue
	 do 124 iz = 1,nz
	   do 124 iy = 1,ny
	     feldx(1,iy,iz) = 0.
	     feldy(1,iy,iz) = 0.
	     feldz(1,iy,iz) = 0.
	     feldx(nx,iy,iz) = 0.
	     feldy(nx,iy,iz) = 0.
	     feldz(nx,iy,iz) = 0.
  124    continue
      return
      end
	subroutine bdcoef
c***********************************************************************
	include 'misflin'
c
	integer  i,j,nx3,ny3,nz3
c ............................................................
c     bdcoef determiines coefficients for the boundary 
c     conditions in the following manner 
c ............................................................
c     f(n) = k*f(n-2) + alfa*f'(n-3)*del
c     mit
c        del = 1/dif(n-1)
c        f'(n-2)  = dif(n-3)*[f(n-2) - f(n-4)]
c   => f(n) = ( k + alfa*dif(n-3)/dif(n-1) )*f(n-2)
c                     - alfa*dif(n-3)/dif(n-1)*f(n-4)
c   und
c      f(1) = ( k + alfa*dif(4)/dif(2) )*f(3)
c                     - alfa*dif(4)/dif(2)*f(5)
c.............................................
c
c     antisymmetry      k = -1  ,  alfa = 0
c     symmetry          k = 1   ,  alfa = 0
c     extrapolation      k = 1   ,  alfa = 1
c
c     b from div b = 0 
c.............................................
c     boundcoeffizienten:
c       1. index:  1=min, 2=max
c       2. index:  1=x, 2=y, 3=z
c
      nx3 = nx-3
      ny3 = ny-3
      nz3 = nz-3
c
c..................
c     x = xmin
c..................
      crho1(1,1) = 1. + arho(1,1)*difx(4)/difx(2)
      crho2(1,1) =  - arho(1,1)*difx(4)/difx(2)
      cu1(1,1)   = 1. + au(1,1)*difx(4)/difx(2)
      cu2(1,1)   =  - au(1,1)*difx(4)/difx(2)
      cres1(1,1)   = 1. + ares(1,1)*difx(4)/difx(2)
      cres2(1,1)   =  - ares(1,1)*difx(4)/difx(2)
      cbx1(1,1)  = kbx(1,1) + abx(1,1)*difx(4)/difx(2)
      cbx2(1,1)  =  - abx(1,1)*difx(4)/difx(2)
      cby1(1,1)  = kby(1,1) + aby(1,1)*difx(4)/difx(2)
      cby2(1,1)  =  - aby(1,1)*difx(4)/difx(2)
      cbz1(1,1)  = kbz(1,1) + abz(1,1)*difx(4)/difx(2)
      cbz2(1,1)  =  - abz(1,1)*difx(4)/difx(2)
      csx1(1,1)  = ksx(1,1) + asx(1,1)*difx(4)/difx(2)
      csx2(1,1)  =  - asx(1,1)*difx(4)/difx(2)
      csy1(1,1)  = ksy(1,1) + asy(1,1)*difx(4)/difx(2)
      csy2(1,1)  =  - asy(1,1)*difx(4)/difx(2)
      csz1(1,1)  = ksz(1,1) + asz(1,1)*difx(4)/difx(2)
      csz2(1,1)  =  - asz(1,1)*difx(4)/difx(2)
c..................
c     x = xmax
c..................
      crho1(2,1) = 1. + arho(2,1)*difx(nx3)/difx(nx1)
      crho2(2,1) =  - arho(2,1)*difx(nx3)/difx(nx1)
      cu1(2,1)   = 1. + au(2,1)*difx(nx3)/difx(nx1)
      cu2(2,1)   =  - au(2,1)*difx(nx3)/difx(nx1)
      cres1(2,1)   = 1. + ares(2,1)*difx(nx3)/difx(nx1)
      cres2(2,1)   =  - ares(2,1)*difx(nx3)/difx(nx1)
      cbx1(2,1)  = kbx(2,1) + abx(2,1)*difx(nx3)/difx(nx1)
      cbx2(2,1)  =  - abx(2,1)*difx(nx3)/difx(nx1)
      cby1(2,1)  = kby(2,1) + aby(2,1)*difx(nx3)/difx(nx1)
      cby2(2,1)  =  - aby(2,1)*difx(nx3)/difx(nx1)
      cbz1(2,1)  = kbz(2,1) + abz(2,1)*difx(nx3)/difx(nx1)
      cbz2(2,1)  =  - abz(2,1)*difx(nx3)/difx(nx1)
      csx1(2,1)  = ksx(2,1) + asx(2,1)*difx(nx3)/difx(nx1)
      csx2(2,1)  =  - asx(2,1)*difx(nx3)/difx(nx1)
      csy1(2,1)  = ksy(2,1) + asy(2,1)*difx(nx3)/difx(nx1)
      csy2(2,1)  =  - asy(2,1)*difx(nx3)/difx(nx1)
      csz1(2,1)  = ksz(2,1) + asz(2,1)*difx(nx3)/difx(nx1)
      csz2(2,1)  =  - asz(2,1)*difx(nx3)/difx(nx1)
c..................
c     y = ymin
c..................
      crho1(1,2) = 1. + arho(1,2)*dify(4)/dify(2)
      crho2(1,2) =  - arho(1,2)*dify(4)/dify(2)
      cu1(1,2)   = 1. + au(1,2)*dify(4)/dify(2)
      cu2(1,2)   =  - au(1,2)*dify(4)/dify(2)
      cres1(1,2)   = 1. + ares(1,2)*dify(4)/dify(2)
      cres2(1,2)   =  - ares(1,2)*dify(4)/dify(2)
      cbx1(1,2)  = kbx(1,2) + abx(1,2)*dify(4)/dify(2)
      cbx2(1,2)  =  - abx(1,2)*dify(4)/dify(2)
      cby1(1,2)  = kby(1,2) + aby(1,2)*dify(4)/dify(2)
      cby2(1,2)  =  - aby(1,2)*dify(4)/dify(2)
      cbz1(1,2)  = kbz(1,2) + abz(1,2)*dify(4)/dify(2)
      cbz2(1,2)  =  - abz(1,2)*dify(4)/dify(2)
      csx1(1,2)  = ksx(1,2) + asx(1,2)*dify(4)/dify(2)
      csx2(1,2)  =  - asx(1,2)*dify(4)/dify(2)
      csy1(1,2)  = ksy(1,2) + asy(1,2)*dify(4)/dify(2)
      csy2(1,2)  =  - asy(1,2)*dify(4)/dify(2)
      csz1(1,2)  = ksz(1,2) + asz(1,2)*dify(4)/dify(2)
      csz2(1,2)  =  - asz(1,2)*dify(4)/dify(2)
c..................
c     y = ymax
c..................
      crho1(2,2) = 1. + arho(2,2)*dify(ny3)/dify(ny1)
      crho2(2,2) =  - arho(2,2)*dify(ny3)/dify(ny1)
      cu1(2,2)   = 1. + au(2,2)*dify(ny3)/dify(ny1)
      cu2(2,2)   =  - au(2,2)*dify(ny3)/dify(ny1)
      cres1(2,2)   = 1. + ares(2,2)*dify(ny3)/dify(ny1)
      cres2(2,2)   =  - ares(2,2)*dify(ny3)/dify(ny1)
      cbx1(2,2)  = kbx(2,2) + abx(2,2)*dify(ny3)/dify(ny1)
      cbx2(2,2)  =  - abx(2,2)*dify(ny3)/dify(ny1)
      cby1(2,2)  = kby(2,2) + aby(2,2)*dify(ny3)/dify(ny1)
      cby2(2,2)  =  - aby(2,2)*dify(ny3)/dify(ny1)
      cbz1(2,2)  = kbz(2,2) + abz(2,2)*dify(ny3)/dify(ny1)
      cbz2(2,2)  =  - abz(2,2)*dify(ny3)/dify(ny1)
      csx1(2,2)  = ksx(2,2) + asx(2,2)*dify(ny3)/dify(ny1)
      csx2(2,2)  =  - asx(2,2)*dify(ny3)/dify(ny1)
      csy1(2,2)  = ksy(2,2) + asy(2,2)*dify(ny3)/dify(ny1)
      csy2(2,2)  =  - asy(2,2)*dify(ny3)/dify(ny1)
      csz1(2,2)  = ksz(2,2) + asz(2,2)*dify(ny3)/dify(ny1)
      csz2(2,2)  =  - asz(2,2)*dify(ny3)/dify(ny1)
c..................
c     z = zmin
c..................
      crho1(1,3) = 1. + arho(1,3)*difz(4)/difz(2)
      crho2(1,3) =  - arho(1,3)*difz(4)/difz(2)
      cu1(1,3)   = 1. + au(1,3)*difz(4)/difz(2)
      cu2(1,3)   =  - au(1,3)*difz(4)/difz(2)
      cres1(1,3)   = 1. + ares(1,3)*difz(4)/difz(2)
      cres2(1,3)   =  - ares(1,3)*difz(4)/difz(2)
      cbx1(1,3)  = kbx(1,3) + abx(1,3)*difz(4)/difz(2)
      cbx2(1,3)  =  - abx(1,3)*difz(4)/difz(2)
      cby1(1,3)  = kby(1,3) + aby(1,3)*difz(4)/difz(2)
      cby2(1,3)  =  - aby(1,3)*difz(4)/difz(2)
      cbz1(1,3)  = kbz(1,3) + abz(1,3)*difz(4)/difz(2)
      cbz2(1,3)  =  - abz(1,3)*difz(4)/difz(2)
      csx1(1,3)  = ksx(1,3) + asx(1,3)*difz(4)/difz(2)
      csx2(1,3)  =  - asx(1,3)*difz(4)/difz(2)
      csy1(1,3)  = ksy(1,3) + asy(1,3)*difz(4)/difz(2)
      csy2(1,3)  =  - asy(1,3)*difz(4)/difz(2)
      csz1(1,3)  = ksz(1,3) + asz(1,3)*difz(4)/difz(2)
      csz2(1,3)  =  - asz(1,3)*difz(4)/difz(2)
c..................
c     z = zmax
c..................
      crho1(2,3) = 1. + arho(2,3)*difz(nz3)/difz(nz1)
      crho2(2,3) =  - arho(2,3)*difz(nz3)/difz(nz1)
      cu1(2,3)   = 1. + au(2,3)*difz(nz3)/difz(nz1)
      cu2(2,3)   =  - au(2,3)*difz(nz3)/difz(nz1)
      cres1(2,3)   = 1. + ares(2,3)*difz(nz3)/difz(nz1)
      cres2(2,3)   =  - ares(2,3)*difz(nz3)/difz(nz1)
      cbx1(2,3)  = kbx(2,3) + abx(2,3)*difz(nz3)/difz(nz1)
      cbx2(2,3)  =  - abx(2,3)*difz(nz3)/difz(nz1)
      cby1(2,3)  = kby(2,3) + aby(2,3)*difz(nz3)/difz(nz1)
      cby2(2,3)  =  - aby(2,3)*difz(nz3)/difz(nz1)
      cbz1(2,3)  = kbz(2,3) + abz(2,3)*difz(nz3)/difz(nz1)
      cbz2(2,3)  =  - abz(2,3)*difz(nz3)/difz(nz1)
      csx1(2,3)  = ksx(2,3) + asx(2,3)*difz(nz3)/difz(nz1)
      csx2(2,3)  =  - asx(2,3)*difz(nz3)/difz(nz1)
      csy1(2,3)  = ksy(2,3) + asy(2,3)*difz(nz3)/difz(nz1)
      csy2(2,3)  =  - asy(2,3)*difz(nz3)/difz(nz1)
      csz1(2,3)  = ksz(2,3) + asz(2,3)*difz(nz3)/difz(nz1)
      csz2(2,3)  =  - asz(2,3)*difz(nz3)/difz(nz1)
c
      do 10 j = 1,3
	if (perio(j)) then
	   lsym(1,j) = .false.
	   lsym(2,j) = .false.
	end if
  10  continue
c
      do 100 i = 1,2
       do 100 j = 1,3
	if ( perio(j) .or. lsym(i,j) ) then
	  crho1(i,j)  = 0.0
	  crho2(i,j)  = 0.0
	  cu1(i,j)    = 0.0
	  cu2(i,j)    = 0.0
	  cres1(i,j)  = 0.0
	  cres2(i,j)  = 0.0
	  cbx1(i,j)   = 0.0
	  cbx2(i,j)   = 0.0
	  cby1(i,j)   = 0.0
	  cby2(i,j)   = 0.0
	  cbz1(i,j)   = 0.0
	  cbz2(i,j)   = 0.0
	  csx1(i,j)   = 0.0
	  csx2(i,j)   = 0.0
	  csy1(i,j)   = 0.0
	  csy2(i,j)   = 0.0
	  csz1(i,j)   = 0.0
	  csz2(i,j)   = 0.0
	end if
  100 continue
      return
      end
	subroutine bound
c***********************************************************************
	include 'misflin'
c
	integer  ix,iz,iy,nx2,nx4,ny2,ny4,nz2,nz4,
     +           nxhalb,nyhalb,nzhalb
c ............................................................
c   computation of boundary values
c       f(i+1) = c1*f(i-1) + c2*f(i-3) + perio*f(n-3)
c                + lsym*f(sym)
c.............................................................
c     (line symmetry or periodicity determ by lsym or perio 
c
      nxhalb = (nx+1)/2
      nyhalb = (ny+1)/2
      nzhalb = (nz+1)/2
      nx2 = nx-2
      nx4 = nx-4
      ny2 = ny-2
      ny4 = ny-4
      nz2 = nz-2
      nz4 = nz-4
c--------------------------------------
c   x = xmin
c       (linesymmetry along z axis)
c--------------------------------------
      if ( perio(1) )  then

       do 110 iz = 2, nz1
       do 110 iy = 2, ny1
	 rho(1,iy,iz) = rho(nx2,iy,iz)
	 u(1,iy,iz)   = u(nx2,iy,iz)
	 res(1,iy,iz) = res(nx2,iy,iz)
	 bx(1,iy,iz)  = bx(nx2,iy,iz)
	 by(1,iy,iz)  = by(nx2,iy,iz)
	 bz(1,iy,iz)  = bz(nx2,iy,iz)
	 sx(1,iy,iz)  = sx(nx2,iy,iz)
	 sy(1,iy,iz)  = sy(nx2,iy,iz)
	 sz(1,iy,iz)  = sz(nx2,iy,iz)
  110  continue

      else if ( lsym(1,1) ) then

       do 120 iz = 2, nz1
       do 120 iy = 2, ny1
	 rho(1,iy,iz) = rho(3,ny+1-iy,iz)
	 u(1,iy,iz)   = u(3,ny+1-iy,iz)
	 res(1,iy,iz) = res(3,ny+1-iy,iz)
	 bx(1,iy,iz)  = bx(3,ny+1-iy,iz)
	 by(1,iy,iz)  = by(3,ny+1-iy,iz)
	 bz(1,iy,iz)  = -bz(3,ny+1-iy,iz)
	 sx(1,iy,iz)  = -sx(3,ny+1-iy,iz)
	 sy(1,iy,iz)  = -sy(3,ny+1-iy,iz)
	 sz(1,iy,iz)  = sz(3,ny+1-iy,iz)
  120  continue
       do 125 iz = 2, nz1
	 bz(2,nyhalb,iz)  = 0.0
	 sx(2,nyhalb,iz)  = 0.0
	 sy(2,nyhalb,iz)  = 0.0
  125  continue

      else

       do 140 iz = 2, nz1
       do 140 iy = 2, ny1
	 rho(1,iy,iz) = crho1(1,1)*rho(3,iy,iz) + crho2(1,1)*rho(5,iy,iz)
	 u(1,iy,iz)   = cu1(1,1)*u(3,iy,iz) + cu2(1,1)*u(5,iy,iz)
	 res(1,iy,iz) = cres1(1,1)*res(3,iy,iz) + cres2(1,1)*res(5,iy,iz)
	 bx(1,iy,iz)  = cbx1(1,1)*bx(3,iy,iz) + cbx2(1,1)*bx(5,iy,iz)
	 by(1,iy,iz)  = cby1(1,1)*by(3,iy,iz) + cby2(1,1)*by(5,iy,iz)
	 bz(1,iy,iz)  = cbz1(1,1)*bz(3,iy,iz) + cbz2(1,1)*bz(5,iy,iz)
	 sx(1,iy,iz)  = csx1(1,1)*sx(3,iy,iz) + csx2(1,1)*sx(5,iy,iz)
	 sy(1,iy,iz)  = csy1(1,1)*sy(3,iy,iz) + csy2(1,1)*sy(5,iy,iz)
	 sz(1,iy,iz)  = csz1(1,1)*sz(3,iy,iz) + csz2(1,1)*sz(5,iy,iz)
  140  continue

       if (ksx(1,1) .eq. -1.) then
	do 141 iz = 2, nz1
	do 141 iy = 2, ny1
  141     sx(2,iy,iz) = 0.
       end if
       if (ksy(1,1) .eq. -1.) then
	do 142 iz = 2, nz1
	do 142 iy = 2, ny1
  142     sy(2,iy,iz) = 0.
       end if
       if (ksz(1,1) .eq. -1.) then
	do 143 iz = 2, nz1
	do 143 iy = 2, ny1
  143     sz(2,iy,iz) = 0.
       end if
       if (kbx(1,1) .eq. -1.) then
	do 146 iz = 2, nz1
	do 146 iy = 2, ny1
  146     bx(2,iy,iz) = 0.
       end if
       if (kby(1,1) .eq. -1.) then
	do 147 iz = 2, nz1
	do 147 iy = 2, ny1
  147     by(2,iy,iz) = 0.
       end if
       if (kbz(1,1) .eq. -1.) then
	do 148 iz = 2, nz1
	do 148 iy = 2, ny1
  148     bz(2,iy,iz) = 0.
       end if


      end if

c--------------------------------------
c   x = xmax
c       (linesymmetry along z axis)
c--------------------------------------
      if ( perio(1) )  then

       do 210 iz = 2, nz1
       do 210 iy = 2, ny1
	rho(nx,iy,iz) = rho(3,iy,iz)
	u(nx,iy,iz)   = u(3,iy,iz)
	res(nx,iy,iz) = res(3,iy,iz)
	bx(nx,iy,iz)  = bx(3,iy,iz)
	by(nx,iy,iz)  = by(3,iy,iz)
	bz(nx,iy,iz)  = bz(3,iy,iz)
	sx(nx,iy,iz)  = sx(3,iy,iz)
	sy(nx,iy,iz)  = sy(3,iy,iz)
	sz(nx,iy,iz)  = sz(3,iy,iz)
  210  continue

      else if ( lsym(2,1) )  then

       do 220 iz = 2, nz1
       do 220 iy = 2, ny1
	rho(nx,iy,iz) =  rho(nx2,ny+1-iy,iz)
	u(nx,iy,iz)   =  u(nx2,ny+1-iy,iz)
	res(nx,iy,iz) =  res(nx2,ny+1-iy,iz)
	bx(nx,iy,iz)  =  bx(nx2,ny+1-iy,iz)
	by(nx,iy,iz)  =  by(nx2,ny+1-iy,iz)
	bz(nx,iy,iz)  = -bz(nx2,ny+1-iy,iz)
	sx(nx,iy,iz)  = -sx(nx2,ny+1-iy,iz)
	sy(nx,iy,iz)  = -sy(nx2,ny+1-iy,iz)
	sz(nx,iy,iz)  =  sz(nx2,ny+1-iy,iz)
  220  continue
       do 225 iz = 2, nz1
	 bz(nx1,nyhalb,iz)  = 0.0
	 sx(nx1,nyhalb,iz)  = 0.0
	 sy(nx1,nyhalb,iz)  = 0.0
  225  continue

      else

       do 240 iz = 2, nz1
       do 240 iy = 2, ny1
	rho(nx,iy,iz)= crho1(2,1)*rho(nx2,iy,iz) 
     +                 + crho2(2,1)*rho(nx4,iy,iz)
	u(nx,iy,iz)  = cu1(2,1)*u(nx2,iy,iz) + cu2(2,1)*u(nx4,iy,iz)
	res(nx,iy,iz)= cres1(2,1)*res(nx2,iy,iz) 
     +                 + cres2(2,1)*res(nx4,iy,iz)
	bx(nx,iy,iz) = cbx1(2,1)*bx(nx2,iy,iz) + cbx2(2,1)*bx(nx4,iy,iz)
	by(nx,iy,iz) = cby1(2,1)*by(nx2,iy,iz) + cby2(2,1)*by(nx4,iy,iz)
	bz(nx,iy,iz) = cbz1(2,1)*bz(nx2,iy,iz) + cbz2(2,1)*bz(nx4,iy,iz)
	sx(nx,iy,iz) = csx1(2,1)*sx(nx2,iy,iz) + csx2(2,1)*sx(nx4,iy,iz)
	sy(nx,iy,iz) = csy1(2,1)*sy(nx2,iy,iz) + csy2(2,1)*sy(nx4,iy,iz)
	sz(nx,iy,iz) = csz1(2,1)*sz(nx2,iy,iz) + csz2(2,1)*sz(nx4,iy,iz)
  240  continue

       if (ksx(2,1) .eq. -1.) then
	do 241 iz = 2, nz1
	do 241 iy = 2, ny1
  241     sx(nx1,iy,iz) = 0.
       end if
       if (ksy(2,1) .eq. -1.) then
	do 242 iz = 2, nz1
	do 242 iy = 2, ny1
  242     sy(nx1,iy,iz) = 0.
       end if
       if (ksz(2,1) .eq. -1.) then
	do 243 iz = 2, nz1
	do 243 iy = 2, ny1
  243     sz(nx1,iy,iz) = 0.
       end if
       if (kbx(2,1) .eq. -1.) then
	do 246 iz = 2, nz1
	do 246 iy = 2, ny1
  246     bx(nx1,iy,iz) = 0.
       end if
       if (kby(2,1) .eq. -1.) then
	do 247 iz = 2, nz1
	do 247 iy = 2, ny1
  247     by(nx1,iy,iz) = 0.
       end if
       if (kbz(2,1) .eq. -1.) then
	do 248 iz = 2, nz1
	do 248 iy = 2, ny1
  248     bz(nx1,iy,iz) = 0.
       end if


      end if

c--------------------------------------
c   y = ymin
c       (linesymmetry along x axis)
c--------------------------------------
      if ( perio(2) )  then

       do 310 iz = 2, nz1
       do 310 ix = 1, nx
	rho(ix,1,iz) = rho(ix,ny2,iz)
	u(ix,1,iz)   = u(ix,ny2,iz)
	res(ix,1,iz) = res(ix,ny2,iz)
	bx(ix,1,iz)  = bx(ix,ny2,iz)
	by(ix,1,iz)  = by(ix,ny2,iz)
	bz(ix,1,iz)  = bz(ix,ny2,iz)
	sx(ix,1,iz)  = sx(ix,ny2,iz)
	sy(ix,1,iz)  = sy(ix,ny2,iz)
	sz(ix,1,iz)  = sz(ix,ny2,iz)
  310  continue

      else if ( lsym(1,2) )  then

       do 320 iz = 2, nz1
       do 320 ix = 1, nx
	rho(ix,1,iz) =  rho(ix,3,nz+1-iz)
	u(ix,1,iz)   =  u(ix,3,nz+1-iz)
	res(ix,1,iz) =  res(ix,3,nz+1-iz)
	bx(ix,1,iz)  = -bx(ix,3,nz+1-iz)
	by(ix,1,iz)  =  by(ix,3,nz+1-iz)
	bz(ix,1,iz)  =  bz(ix,3,nz+1-iz)
	sx(ix,1,iz)  =  sx(ix,3,nz+1-iz)
	sy(ix,1,iz)  = -sy(ix,3,nz+1-iz)
	sz(ix,1,iz)  = -sz(ix,3,nz+1-iz)
  320 continue
       do 325 ix = 1, nx
	bx(ix,2,nzhalb)  = 0.0
	sy(ix,2,nzhalb)  = 0.0
	sz(ix,2,nzhalb)  = 0.0
  325  continue

      else

       do 340 iz = 2, nz1
       do 340 ix = 1, nx
	rho(ix,1,iz) = crho1(1,2)*rho(ix,3,iz)
     +                + crho2(1,2)*rho(ix,5,iz)
	u(ix,1,iz)   = cu1(1,2)*u(ix,3,iz) + cu2(1,2)*u(ix,5,iz)
	res(ix,1,iz)   = cres1(1,2)*res(ix,3,iz)
     +                + cres2(1,2)*res(ix,5,iz)
	bx(ix,1,iz)  = cbx1(1,2)*bx(ix,3,iz)+cbx2(1,2)*bx(ix,5,iz)
	by(ix,1,iz)  = cby1(1,2)*by(ix,3,iz)+cby2(1,2)*by(ix,5,iz)
	bz(ix,1,iz)  = cbz1(1,2)*bz(ix,3,iz)+cbz2(1,2)*bz(ix,5,iz)
	sx(ix,1,iz)  = csx1(1,2)*sx(ix,3,iz)+csx2(1,2)*sx(ix,5,iz)
	sy(ix,1,iz)  = csy1(1,2)*sy(ix,3,iz)+csy2(1,2)*sy(ix,5,iz)
	sz(ix,1,iz)  = csz1(1,2)*sz(ix,3,iz)+csz2(1,2)*sz(ix,5,iz)
  340  continue

       if (ksx(1,2) .eq. -1.) then
	do 341 iz = 2, nz1
	do 341 ix = 1, nx
  341     sx(ix,2,iz) = 0.
       end if
       if (ksy(1,2) .eq. -1.) then
	do 342 iz = 2, nz1
	do 342 ix = 1, nx
  342     sy(ix,2,iz) = 0.
       end if
       if (ksz(1,2) .eq. -1.) then
	do 343 iz = 2, nz1
	do 343 ix = 1, nx
  343     sz(ix,2,iz) = 0.
       end if
       if (kbx(1,2) .eq. -1.) then
	do 346 iz = 2, nz1
	do 346 ix = 1, nx
  346     bx(ix,2,iz) = 0.
       end if
       if (kby(1,2) .eq. -1.) then
	do 347 iz = 2, nz1
	do 347 ix = 1, nx
  347     by(ix,2,iz) = 0.
       end if
       if (kbz(1,2) .eq. -1.) then
	do 348 iz = 2, nz1
	do 348 ix = 1, nx
  348     bz(ix,2,iz) = 0.
       end if


      end if

c--------------------------------------
c   y = ymax
c       (linesymmetry along x axis)
c--------------------------------------
      if ( perio(2) )  then

       do 410 iz = 2, nz1
       do 410 ix = 1, nx
	rho(ix,ny,iz) = rho(ix,3,iz)
	u(ix,ny,iz)   = u(ix,3,iz)
	res(ix,ny,iz) = res(ix,3,iz)
	bx(ix,ny,iz)  = bx(ix,3,iz)
	by(ix,ny,iz)  = by(ix,3,iz)
	bz(ix,ny,iz)  = bz(ix,3,iz)
	sx(ix,ny,iz)  = sx(ix,3,iz)
	sy(ix,ny,iz)  = sy(ix,3,iz)
	sz(ix,ny,iz)  = sz(ix,3,iz)
  410  continue

      else if ( lsym(2,2) )  then

       do 420 iz = 2, nz1
       do 420 ix = 1, nx
	rho(ix,ny,iz) =  rho(ix,ny2,nz+1-iz)
	u(ix,ny,iz)   =  u(ix,ny2,nz+1-iz)
	res(ix,ny,iz) =  res(ix,ny2,nz+1-iz)
	bx(ix,ny,iz)  = -bx(ix,ny2,nz+1-iz)
	by(ix,ny,iz)  =  by(ix,ny2,nz+1-iz)
	bz(ix,ny,iz)  =  bz(ix,ny2,nz+1-iz)
	sx(ix,ny,iz)  =  sx(ix,ny2,nz+1-iz)
	sy(ix,ny,iz)  = -sy(ix,ny2,nz+1-iz)
	sz(ix,ny,iz)  = -sz(ix,ny2,nz+1-iz)
  420  continue
       do 425 ix = 1, nx
	bx(ix,ny1,nzhalb)  = 0.0
	sy(ix,ny1,nzhalb)  = 0.0
	sz(ix,ny1,nzhalb)  = 0.0
  425  continue

      else

       do 440 iz = 2, nz1
       do 440 ix = 1, nx
	rho(ix,ny,iz) = crho1(2,2)*rho(ix,ny2,iz)
     +                + crho2(2,2)*rho(ix,ny4,iz)
	u(ix,ny,iz)   = cu1(2,2)*u(ix,ny2,iz)+cu2(2,2)*u(ix,ny4,iz)
	res(ix,ny,iz) = cres1(2,2)*res(ix,ny2,iz)
     +                + cres2(2,2)*res(ix,ny4,iz)
	bx(ix,ny,iz)  = cbx1(2,2)*bx(ix,ny2,iz)+cbx2(2,2)*bx(ix,ny4,iz)
	by(ix,ny,iz)  = cby1(2,2)*by(ix,ny2,iz)+cby2(2,2)*by(ix,ny4,iz)
	bz(ix,ny,iz)  = cbz1(2,2)*bz(ix,ny2,iz)+cbz2(2,2)*bz(ix,ny4,iz)
	sx(ix,ny,iz)  = csx1(2,2)*sx(ix,ny2,iz)+csx2(2,2)*sx(ix,ny4,iz)
	sy(ix,ny,iz)  = csy1(2,2)*sy(ix,ny2,iz)+csy2(2,2)*sy(ix,ny4,iz)
	sz(ix,ny,iz)  = csz1(2,2)*sz(ix,ny2,iz)+csz2(2,2)*sz(ix,ny4,iz)
  440  continue

       if (ksx(2,2) .eq. -1.) then
	do 441 iz = 2, nz1
	do 441 ix = 1, nx
  441     sx(ix,ny1,iz) = 0.
       end if
       if (ksy(2,2) .eq. -1.) then
	do 442 iz = 2, nz1
	do 442 ix = 1, nx
  442     sy(ix,ny1,iz) = 0.
       end if
       if (ksz(2,2) .eq. -1.) then
	do 443 iz = 2, nz1
	do 443 ix = 1, nx
  443     sz(ix,ny1,iz) = 0.
       end if
       if (kbx(2,2) .eq. -1.) then
	do 446 iz = 2, nz1
	do 446 ix = 1, nx
  446     bx(ix,ny1,iz) = 0.
       end if
       if (kby(2,2) .eq. -1.) then
	do 447 iz = 2, nz1
	do 447 ix = 1, nx
  447     by(ix,ny1,iz) = 0.
       end if
       if (kbz(2,2) .eq. -1.) then
	do 448 iz = 2, nz1
	do 448 ix = 1, nx
  448     bz(ix,ny1,iz) = 0.
       end if


      end if

c--------------------------------------
c   z = zmin
c--------------------------------------
      if ( perio(3) )  then

       do 510 iy = 1, ny
       do 510 ix = 1, nx
	rho(ix,iy,1) = rho(ix,iy,nz2)
	u(ix,iy,1)   = u(ix,iy,nz2)
	res(ix,iy,1) = res(ix,iy,nz2)
	bx(ix,iy,1)  = bx(ix,iy,nz2)
	by(ix,iy,1)  = by(ix,iy,nz2)
	bz(ix,iy,1)  = bz(ix,iy,nz2)
	sx(ix,iy,1)  = sx(ix,iy,nz2)
	sy(ix,iy,1)  = sy(ix,iy,nz2)
	sz(ix,iy,1)  = sz(ix,iy,nz2)
  510  continue

      else if ( lsym(1,3) )  then

c for magnetopause computation linesymmetry along x
       do 520 iy = 1, ny
       do 520 ix = 1, nx
	rho(ix,iy,1) =  rho(ix,ny+1-iy,3)
	u(ix,iy,1)   =  u(ix,ny+1-iy,3)
	res(ix,iy,1) =  res(ix,ny+1-iy,3)
	bx(ix,iy,1)  = -bx(ix,ny+1-iy,3)
	by(ix,iy,1)  =  by(ix,ny+1-iy,3)
	bz(ix,iy,1)  =  bz(ix,ny+1-iy,3)
	sx(ix,iy,1)  =  sx(ix,ny+1-iy,3)
	sy(ix,iy,1)  = -sy(ix,ny+1-iy,3)
	sz(ix,iy,1)  = -sz(ix,ny+1-iy,3)
  520  continue
       do 525 ix = 1,nx
	bx(ix,nyhalb,2)  = 0.0
	sy(ix,nyhalb,2)  = 0.0
	sz(ix,nyhalb,2)  = 0.0
  525  continue

c for arcs ( y(pro)->z, z(pro)->-y ) linesym along y(pro):
c       do 520 iy = 1, ny
c       do 520 ix = 1, nx
c	rho(ix,iy,1) =  rho(nx+1-ix,iy,3)
c	u(ix,iy,1)   =  u(nx+1-ix,iy,3)
c	res(ix,iy,1) =  res(nx+1-ix,iy,3)
c	bx(ix,iy,1)  = -bx(nx+1-ix,iy,3)
c	by(ix,iy,1)  =  by(nx+1-ix,iy,3)
c	bz(ix,iy,1)  = -bz(nx+1-ix,iy,3)
c	sx(ix,iy,1)  = -sx(nx+1-ix,iy,3)
c	sy(ix,iy,1)  =  sy(nx+1-ix,iy,3)
c	sz(ix,iy,1)  = -sz(nx+1-ix,iy,3)
c  520  continue
c       do 525 iy = 1,ny
c	bx(nxhalb,iy,2)  = 0.0
c	bz(nxhalb,iy,2)  = 0.0
c	sx(nxhalb,iy,2)  = 0.0
c	sz(nxhalb,iy,2)  = 0.0
c  525  continue

      else

       do 540 iy = 1, ny
       do 540 ix = 1, nx
	rho(ix,iy,1)= crho1(1,3)*rho(ix,iy,3) + crho2(1,3)*rho(ix,iy,5)
	u(ix,iy,1)  = cu1(1,3)*u(ix,iy,3) + cu2(1,3)*u(ix,iy,5)
	res(ix,iy,1)= cres1(1,3)*res(ix,iy,3) + cres2(1,3)*res(ix,iy,5)
	bx(ix,iy,1) = cbx1(1,3)*bx(ix,iy,3) + cbx2(1,3)*bx(ix,iy,5)
	by(ix,iy,1) = cby1(1,3)*by(ix,iy,3) + cby2(1,3)*by(ix,iy,5)
	bz(ix,iy,1) = cbz1(1,3)*bz(ix,iy,3) + cbz2(1,3)*bz(ix,iy,5)
	sx(ix,iy,1) = csx1(1,3)*sx(ix,iy,3) + csx2(1,3)*sx(ix,iy,5)
	sy(ix,iy,1) = csy1(1,3)*sy(ix,iy,3) + csy2(1,3)*sy(ix,iy,5)
	sz(ix,iy,1) = csz1(1,3)*sz(ix,iy,3) + csz2(1,3)*sz(ix,iy,5)
  540  continue

       if (ksx(1,3) .eq. -1.) then
	do 541 iy = 1, ny
	do 541 ix = 1, nx
  541     sx(ix,iy,2) = 0.
       end if
       if (ksy(1,3) .eq. -1.) then
	do 542 iy = 1, ny
	do 542 ix = 1, nx
  542     sy(ix,iy,2) = 0.
       end if
       if (ksz(1,3) .eq. -1.) then
	do 543 iy = 1, ny
	do 543 ix = 1, nx
  543     sz(ix,iy,2) = 0.
       end if
       if (kbx(1,3) .eq. -1.) then
	do 546 iy = 1, ny
	do 546 ix = 1, nx
  546     bx(ix,iy,2) = 0.
       end if
       if (kby(1,3) .eq. -1.) then
	do 547 iy = 1, ny
	do 547 ix = 1, nx
  547     by(ix,iy,2) = 0.
       end if
       if (kbz(1,3) .eq. -1.) then
	do 548 iy = 1, ny
	do 548 ix = 1, nx
  548     bz(ix,iy,2) = 0.
       end if


      end if

c--------------------------------------
c   z = zmax
c--------------------------------------
      if ( perio(3) )  then

       do 610 iy = 1, ny
       do 610 ix = 1, nx
	rho(ix,iy,nz) = rho(ix,iy,3)
	u(ix,iy,nz)   = u(ix,iy,3)
	res(ix,iy,nz) = res(ix,iy,3)
	bx(ix,iy,nz)  = bx(ix,iy,3)
	by(ix,iy,nz)  = by(ix,iy,3)
	bz(ix,iy,nz)  = bz(ix,iy,3)
	sx(ix,iy,nz)  = sx(ix,iy,3)
	sy(ix,iy,nz)  = sy(ix,iy,3)
	sz(ix,iy,nz)  = sz(ix,iy,3)
  610  continue

      else if ( lsym(2,3) )  then

c for magnetopause computation linesymmetry along x
       do 620 iy = 1, ny
       do 620 ix = 1, nx
	rho(ix,iy,nz) =  rho(ix,ny+1-iy,nz2)
	u(ix,iy,nz)   =  u(ix,ny+1-iy,nz2)
	res(ix,iy,nz) =  res(ix,ny+1-iy,nz2)
	bx(ix,iy,nz)  = -bx(ix,ny+1-iy,nz2)
	by(ix,iy,nz)  =  by(ix,ny+1-iy,nz2)
	bz(ix,iy,nz)  =  bz(ix,ny+1-iy,nz2)
	sx(ix,iy,nz)  =  sx(ix,ny+1-iy,nz2)
	sy(ix,iy,nz)  = -sy(ix,ny+1-iy,nz2)
	sz(ix,iy,nz)  = -sz(ix,ny+1-iy,nz2)
  620  continue
       do 625 ix = 1,nx
	bx(ix,nyhalb,nz1)  = 0.0
	sy(ix,nyhalb,nz1)  = 0.0
	sz(ix,nyhalb,nz1)  = 0.0
  625  continue

c       do 620 iy = 1, ny
c       do 620 ix = 1, nx
c	rho(ix,iy,nz) =  rho(nx+1-ix,iy,nz2)
c	u(ix,iy,nz)   =  u(nx+1-ix,iy,nz2)
c	res(ix,iy,nz) =  res(nx+1-ix,iy,nz2)
c	bx(ix,iy,nz)  = -bx(nx+1-ix,iy,nz2)
c	by(ix,iy,nz)  =  by(nx+1-ix,iy,nz2)
c	bz(ix,iy,nz)  = -bz(nx+1-ix,iy,nz2)
c	sx(ix,iy,nz)  = -sx(nx+1-ix,iy,nz2)
c	sy(ix,iy,nz)  =  sy(nx+1-ix,iy,nz2)
c	sz(ix,iy,nz)  = -sz(nx+1-ix,iy,nz2)
c  620  continue
c       do 625 iy = 1,ny
c	bx(nxhalb,iy,nz1)  = 0.0
c	bz(nxhalb,iy,nz1)  = 0.0
c	sx(nxhalb,iy,nz1)  = 0.0
c	sz(nxhalb,iy,nz1)  = 0.0
c  625  continue

      else

       do 640 iy = 1, ny
       do 640 ix = 1, nx
	rho(ix,iy,nz)= crho1(2,3)*rho(ix,iy,nz2) 
     +                 + crho2(2,3)*rho(ix,iy,nz4)
	u(ix,iy,nz)  = cu1(2,3)*u(ix,iy,nz2) + cu2(2,3)*u(ix,iy,nz4)
	res(ix,iy,nz)= cres1(2,3)*res(ix,iy,nz2) 
     +                 + cres2(2,3)*res(ix,iy,nz4)
	bx(ix,iy,nz) = cbx1(2,3)*bx(ix,iy,nz2) + cbx2(2,3)*bx(ix,iy,nz4)
	by(ix,iy,nz) = cby1(2,3)*by(ix,iy,nz2) + cby2(2,3)*by(ix,iy,nz4)
	bz(ix,iy,nz) = cbz1(2,3)*bz(ix,iy,nz2) + cbz2(2,3)*bz(ix,iy,nz4)
	sx(ix,iy,nz) = csx1(2,3)*sx(ix,iy,nz2) + csx2(2,3)*sx(ix,iy,nz4)
	sy(ix,iy,nz) = csy1(2,3)*sy(ix,iy,nz2) + csy2(2,3)*sy(ix,iy,nz4)
	sz(ix,iy,nz) = csz1(2,3)*sz(ix,iy,nz2) + csz2(2,3)*sz(ix,iy,nz4)
  640  continue

       if (ksx(2,3) .eq. -1.) then
	do 641 iy = 1, ny
	do 641 ix = 1, nx
  641     sx(ix,iy,nz1) = 0.
       end if
       if (ksy(2,3) .eq. -1.) then
	do 642 iy = 1, ny
	do 642 ix = 1, nx
  642     sy(ix,iy,nz1) = 0.
       end if
       if (ksz(2,3) .eq. -1.) then
	do 643 iy = 1, ny
	do 643 ix = 1, nx
  643     sz(ix,iy,nz1) = 0.
       end if
       if (kbx(2,3) .eq. -1.) then
	do 646 iy = 1, ny
	do 646 ix = 1, nx
  646     bx(ix,iy,nz1) = 0.
       end if
       if (kby(2,3) .eq. -1.) then
	do 647 iy = 1, ny
	do 647 ix = 1, nx
  647     by(ix,iy,nz1) = 0.
       end if
       if (kbz(2,3) .eq. -1.) then
	do 648 iy = 1, ny
	do 648 ix = 1, nx
  648     bz(ix,iy,nz1) = 0.
       end if


      end if
c------------------------------------------------
c  div b=0 ! attention at edges
c------------------------------------------------
c
c     div b for bx 
c
c      a) bound cond at xmin,xmax
c      b) edges x,y = const
c      c) corners at zmin,zmax
c------------------------------------------------
c
c   x = xmin bzw xmax
c................................................
      if ( .not.(perio(1) .or. lsym(1,1)) )  then
       do 810 iz = 2, nz1
       do 810 iy = 2, ny1
	  bx(1,iy,iz)  = bx(3,iy,iz) + 1./difx(2) * (
     +                + dify(iy)*(by(3,iy+1,iz)-by(3,iy-1,iz))
     +                + difz(iz)*(bz(3,iy,iz+1)-bz(3,iy,iz-1)) )
  810  continue
      end if

      if ( .not.(perio(1) .or. lsym(2,1)) )  then
       do 820 iz = 2, nz1
       do 820 iy = 2, ny1
	  bx(nx,iy,iz) = bx(nx2,iy,iz) - 1./difx(nx1) * (
     +                + dify(iy)*(by(nx1,iy+1,iz)-by(nx1,iy-1,iz))
     +                + difz(iz)*(bz(nx1,iy,iz+1)-bz(nx1,iy,iz-1)) )
  820  continue
      end if
c
c   y = ymin bzw ymax
c................................................
      if ( .not.(perio(2) .or. lsym(1,2)) )  then
       do 830 iz = 2, nz1
       do 830 ix = 2, nx1
	  by(ix,1,iz) = by(ix,3,iz) + 1./dify(2) * (
     +                + difz(iz)*(bz(ix,2,iz+1)-bz(ix,2,iz-1))
     +                + difx(ix)*(bx(ix+1,2,iz)-bx(ix-1,2,iz)) )
  830  continue
      end if

      if ( .not.(perio(2) .or. lsym(2,2)) )  then
       do 840 iz = 2, nz1
       do 840 ix = 2, nx1
	  by(ix,ny,iz) = by(ix,ny2,iz) - 1./dify(ny1) * (
     +                + difz(iz)*(bz(ix,ny1,iz+1)-bz(ix,ny1,iz-1))
     +                + difx(ix)*(bx(ix+1,ny1,iz)-bx(ix-1,ny1,iz)) )
  840  continue
      end if
c
c    z = zmin bzw zmax
c................................................
      if ( .not.(perio(3) .or. lsym(1,3)) )  then
       do 860 iy = 2, ny1
       do 860 ix = 2, nx1
	 bz(ix,iy,1)  = bz(ix,iy,3) + 1./difz(2) * (
     +               + difx(ix)*(bx(ix+1,iy,2)-bx(ix-1,iy,2))
     +               + dify(iy)*(by(ix,iy+1,2)-by(ix,iy-1,2)) )
  860  continue
      end if

      if ( .not.(perio(3) .or. lsym(2,3)) )  then
       do 862 iy = 2, ny1
       do 862 ix = 2, nx1
	 bz(ix,iy,nz)  = bz(ix,iy,nz2) - 1./difz(nz1) * (
     +                + difx(ix)*(bx(ix+1,iy,nz1)-bx(ix-1,iy,nz1))
     +                + dify(iy)*(by(ix,iy+1,nz1)-by(ix,iy-1,nz1)) )
  862  continue
      end if
c
      return
      end
	subroutine resbd
c***********************************************************************
	include 'misflin'
c
	integer  ix,iy,iz,nx2,nx4,ny2,ny4,nz2,nz4,
     +           nxhalb,nyhalb,nzhalb
c ............................................................
c   computation of boundary values for the resist
c       f(i+1) = c1*f(i-1) + c2*f(i-3) + perio*f(n-3)
c                + lsym*f(sym)
c.............................................................
c     describtion in subroutine rdkoef
c     (line symmetry parameter lsym)
c
      nxhalb = (nx+1)/2
      nyhalb = (ny+1)/2
      nzhalb = (nz+1)/2
      nx2 = nx-2
      nx4 = nx-4
      ny2 = ny-2
      ny4 = ny-4
      nz2 = nz-2
      nz4 = nz-4
c--------------------------------------
c   x = xmin
c       (linesymmetry along z axis)
c--------------------------------------
      if ( perio(1) )  then
       do 110 iz = 2, nz1
       do 110 iy = 2, ny1
	 res(1,iy,iz) = res(3,iy,iz)
  110  continue
      else if ( lsym(1,1) ) then
       do 120 iz = 2, nz1
       do 120 iy = 2, ny1
	 res(1,iy,iz) = res(3,ny+1-iy,iz)
  120  continue
      else
       do 140 iz = 2, nz1
       do 140 iy = 2, ny1
	 res(1,iy,iz) = cres1(1,1)*res(3,iy,iz) + cres2(1,1)*res(5,iy,iz)
  140  continue
      end if
c--------------------------------------
c   x = xmax
c       (linesymmetry along z axis)
c--------------------------------------
      if ( perio(1) )  then
       do 210 iz = 2, nz1
       do 210 iy = 2, ny1
	res(nx,iy,iz) = res(3,iy,iz)
  210  continue
      else if ( lsym(2,1) )  then
       do 220 iz = 2, nz1
       do 220 iy = 2, ny1
	res(nx,iy,iz) =  res(nx2,ny+1-iy,iz)
  220  continue
      else
       do 240 iz = 2, nz1
       do 240 iy = 2, ny1
	res(nx,iy,iz) = cres1(2,1)*res(nx2,iy,iz) 
     +                  + cres2(2,1)*res(nx4,iy,iz)
  240  continue
      end if
c--------------------------------------
c   y = ymin
c       (linesymmetry along x axis)
c--------------------------------------
      if ( perio(2) )  then
       do 310 iz = 2, nz1
       do 310 ix = 1, nx
	res(ix,1,iz) = res(ix,ny2,iz)
  310  continue
      else if ( lsym(1,2) )  then
       do 320 iz = 2, nz1
       do 320 ix = 1, nx
	res(ix,1,iz) =  res(ix,3,nz+1-iz)
  320 continue
      else
       do 340 iz = 2, nz1
       do 340 ix = 1, nx
	res(ix,1,iz)   = cres1(1,2)*res(ix,3,iz)
     +                + cres2(1,2)*res(ix,5,iz)
  340  continue
      end if
c--------------------------------------
c   y = ymax
c       (linesymmetry along x axis)
c--------------------------------------
      if ( perio(2) )  then
       do 410 iz = 2, nz1
       do 410 ix = 1, nx
	res(ix,ny,iz) = res(ix,3,iz)
  410  continue
      else if ( lsym(2,2) )  then
       do 420 iz = 2, nz1
       do 420 ix = 1, nx
	res(ix,ny,iz) =  res(ix,ny2,nz+1-iz)
  420  continue
      else
       do 440 iz = 2, nz1
       do 440 ix = 1, nx
	res(ix,ny,iz) = cres1(2,2)*res(ix,ny2,iz)
     +                + cres2(2,2)*res(ix,ny4,iz)
  440  continue
      end if
c--------------------------------------
c   z = zmin
c--------------------------------------
      if ( perio(3) )  then
       do 510 iy = 1, ny
       do 510 ix = 1, nx
	res(ix,iy,1) = res(ix,iy,nz2)
  510  continue
      else if ( lsym(1,3) )  then
       do 520 iy = 1, ny
       do 520 ix = 1, nx
	res(ix,iy,1) =  res(ix,ny+1-iy,3)
c	res(ix,iy,1) =  res(nx+1-ix,iy,3)
  520  continue
      else
       do 540 iy = 1, ny
       do 540 ix = 1, nx
	res(ix,iy,1) = cres1(1,3)*res(ix,iy,3) + cres2(1,3)*res(ix,iy,5)
  540  continue
      end if
c--------------------------------------
c   z = zmax
c--------------------------------------
      if ( perio(3) )  then
       do 610 iy = 1, ny
       do 610 ix = 1, nx
	res(ix,iy,nz) = res(ix,iy,3)
  610  continue
      else if ( lsym(2,3) )  then
       do 620 iy = 1, ny
       do 620 ix = 1, nx
	res(ix,iy,nz) =  res(ix,ny+1-iy,nz2)
c	res(ix,iy,nz) =  res(nx+1-ix,iy,nz2)
  620  continue
      else
       do 640 iy = 1, ny
       do 640 ix = 1, nx
	res(ix,iy,nz) = cres1(2,3)*res(ix,iy,nz2) 
     +                  + cres2(2,3)*res(ix,iy,nz4)
  640  continue
      end if
c
      return
      end
	subroutine termin
c***********************************************************************
	include 'misflin'
c
	integer  ix,iy,iz,ksum,tstabx(nx),tstaby(ny),tstabz(nz),
     +           ixphase,iyphase,izphase,ixgrid,iygrid,izgrid,
     +           indx(nx),indy(ny),indz(nz),is,iymax(nx),izmax(nx),
     +           nox,noy,noz
	real     absv2,vphase,vpgrid,dgrid,vp(nx)
c ........................................................
c  termination of program for rho,u .lt. 0 oder v .gt. 4
c ........................................................
      do 20 ix = 1,nx
         tstabx(ix) = 0
   20 continue
      do 22 iy = 1,ny
         tstaby(iy) = 0
   22 continue
      do 24 iz = 1,nz
         tstabz(iz) = 0
   24 continue
      do 30 iz = 1,nz
      do 30 iy = 1,ny
      do 30 ix = 1,nx
       if ( rho(ix,iy,iz) .le. 0.1 .or.
     +      u(ix,iy,iz) .le. 0.1 .or.
     +      ( sx(ix,iy,iz)*sx(ix,iy,iz) + sy(ix,iy,iz)*sy(ix,iy,iz)
     +       + sz(ix,iy,iz)*sz(ix,iy,iz) )/rho(ix,iy,iz)*rho(ix,iy,iz)
     +                   .ge. 16.0 )  then 
          tstabx(ix) = 1
          tstaby(iy) = 1
          tstabz(iz) = 1
        endif
   30 continue
      ksum = 0
      do 50 ix = 1,nx
        ksum = ksum + tstabx(ix) 
   50 continue
      if (ksum .ge. 1) then
        ferror = .true.
        write(26,38) ksum,istep
        write(*,38) ksum,istep
        is = 0
        do 60 ix = 1,nx
           is = is + tstabx(ix)
           if (tstabx(ix).eq.1) indx(is) = ix 
   60   continue
        nox = is
        is = 0
        do 70 iy = 1,ny
           is = is + tstaby(iy)
           if (tstaby(iy).eq.1) indy(is) = iy 
   70   continue
        noy = is
        is = 0
        do 80 iz = 1,nz
           is = is + tstabz(iz)
           if (tstabz(iz).eq.1) then
              indz(is) = iz 
              write(*,*) 'step 2a, is, indz=', is,indz(is) 
           endif
   80   continue
        noz = is
   
       	do 100 iz = 1,noz
       	do 100 iy = 1,noy
        do 100 ix = 1,nox
	  if ( rho(indx(ix),indy(iy),indz(iz)) .le. 0.1 ) then
	    write(26,39) indx(ix),indy(iy),indz(iz),
     +                   rho(indx(ix),indy(iy),indz(iz)),time
	    rho(indx(ix),indy(iy),indz(iz)) = rhoprof(indx(ix))
          endif
	  if ( u(indx(ix),indy(iy),indz(iz)) .le. 0.1 ) then
	    write(26,40) indx(ix),indy(iy),indz(iz),
     +                   u(indx(ix),indy(iy),indz(iz)),time
	    u(indx(ix),indy(iy),indz(iz)) = uprof(indx(ix))
          endif
          absv2 = ( sx(indx(ix),indy(iy),indz(iz))
     +              *sx(indx(ix),indy(iy),indz(iz))
     +             +sy(indx(ix),indy(iy),indz(iz))
     +              *sy(indx(ix),indy(iy),indz(iz))
     +             +sz(indx(ix),indy(iy),indz(iz))
     +              *sz(indx(ix),indy(iy),indz(iz)) )
     +               /rho(indx(ix),indy(iy),indz(iz))
     +               /rho(indx(ix),indy(iy),indz(iz))
	  if ( absv2.gt.25.0) then
     	     write(26,41) indx(ix),indy(iy),indz(iz),absv2,time
          end if
  100   continue
        do 300 iz = 1,nz
        do 300 iy = 1,ny
        do 300 ix = 1,nx
          help(ix,iy,iz) = sqrt ( ( bx(ix,iy,iz)*bx(ix,iy,iz)
     +         + by(ix,iy,iz)*by(ix,iy,iz) + bz(ix,iy,iz)*bz(ix,iy,iz)
     +         + gamma*u(ix,iy,iz)**gamma ) / rho(ix,iy,iz) )
          hilf(ix,iy,iz) = max(difx(ix),dify(iy),difz(iz))
  300   continue
        vphase = 0.
        vpgrid = 0.
        do 350 ix = 1,nx
           vp(ix) = help(ix,1,1)
           iymax(ix) = 2
           izmax(ix) = 2
  350   continue
        do 360 iz = 1,nz
        do 360 iy = 2,ny 
        do 360 ix = 1,nx 
          if (help(ix,iy,iz).gt.vp(ix)) then
            vp(ix) = help(ix,iy,iz)
            iymax(ix) = iy
            izmax(ix) = iz
          endif
  360   continue
        do 370 ix = 1,nx 
          if (vp(ix).gt.vphase) then
           vphase = vp(ix)
           ixphase=ix
           iyphase=iymax(ix)
           izphase=izmax(ix) 
          endif    
  370   continue
        do 380 iz = 1,nz
        do 380 iy = 1,ny
        do 380 ix = 1,nx
           help(ix,iy,iz) = help(ix,iy,iz)*hilf(ix,iy,iz)
  380   continue
        do 390 ix = 1,nx
           vp(ix) = help(ix,1,1)
           iymax(ix) = 2
           izmax(ix) = 2
  390   continue
        do 400 iz = 1,nz
        do 400 iy = 2,ny 
        do 400 ix = 1,nx 
          if (help(ix,iy,iz).gt.vp(ix)) then
            vp(ix) = help(ix,iy,iz)
            iymax(ix) = iy
            izmax(ix) = iz
          endif
  400   continue
        do 410 ix = 1,nx 
          if (vp(ix).gt.vpgrid) then
           vpgrid = vp(ix)
           ixgrid=ix
           iygrid=iymax(ix)
           izgrid=izmax(ix)
          endif    
  410   continue
        vpgrid=2.*vpgrid
        write(26,42) ixphase,iyphase,izphase,vphase,
     +                ixgrid,iygrid,izgrid,vpgrid,dt 
        if (dt .gt. 0.9/vpgrid) then
           newdt = .true.
           write(26,43) time, 0.9*dt 
        end if
      end if
   38 format(/,'No of errrors is greater than ',i3,'  at istep:',i5)
   39 format(' rho(',i3,',',i3,',',i3,') =',f12.7,'   time =',f9.4)
   40 format(' u(',i3,',',i3,',',i3,')   =',f12.7,'   time =',f9.4)
   41 format(' abs v**2(',i3,',',i3,',',i3,') =',f12.7,
     +       '   time =',f9.4)
   42 format(' vphase(',i3,',',i3,',',i3,') =',f8.3,
     +       '   vpgrid(',i3,',',i3,',',i3,') =',f7.2,
     +       ' present dt =',f6.4)
   43 format(' time =',f10.5,'   newdt =', f7.4)
c
      return
      end
	subroutine diag0
c***********************************************************************
	include 'misflin'
c
        integer  iz,ix,iy
c.....................................................
c   calculate energy density
c.....................................................
      do 40 iz = 1,nz
      do 40 iy = 1,ny
      do 40 ix = 1,nx
        help(ix,iy,iz)  = u(ix,iy,iz)**gamma
	heb(ix,iy,iz)   = 0.5*( bx(ix,iy,iz)*bx(ix,iy,iz)
     +                         + by(ix,iy,iz)*by(ix,iy,iz)
     +                          + bz(ix,iy,iz)*bz(ix,iy,iz) )
	hev(ix,iy,iz)   = 0.5*( sx(ix,iy,iz)*sx(ix,iy,iz)
     +                    + sy(ix,iy,iz)*sy(ix,iy,iz)
     +                     + sz(ix,iy,iz)*sz(ix,iy,iz) )/rho(ix,iy,iz) 
   40 continue
      return
      end
 	subroutine diag1
c***********************************************************************
	include 'misflin'
c
        integer  iz,ix,iy,nxe,icomp,iymax(nx),izmax(nx),
     +           iymin(nx),izmin(nx)
        real     bmax(4,4), bmin(4,4), vmax(4,4), vmin(4,4),
     +           jmax(4,4), jmin(4,4), emax(4,4), emin(4,4),
     +           jemin(2,4),jemax(2,4),prho(4,4),
     +           fmax(nx),fmin(nx)

c.....................................................
c   mass und energy of system
c.....................................................
      nxe = (nx+1)/2
      do 20 iz = 1,nz
      do 20 iy = 1,ny
      do 20 ix = 1,nx
	hilf(ix,iy,iz) = 0.125/difx(ix)/dify(iy)/difz(iz)
	if (ix .eq. 2 .or. ix .eq. nx1) 
     +          hilf(ix,iy,iz) = 0.5*hilf(ix,iy,iz) 
	if (iy .eq. 2 .or. iy .eq. ny1) 
     +          hilf(ix,iy,iz) = 0.5*hilf(ix,iy,iz) 
	if (iz .eq. 2 .or. iz .eq. nz1) 
     +          hilf(ix,iy,iz) = 0.5*hilf(ix,iy,iz) 
   20 continue   

      do 50 iz = 2,nz-1
      do 50 iy = 2,ny-1
      do 50 ix = 2,nx-1
        feldx(ix,iy,iz) = dify(iy)*(bz(ix,iy+1,iz)-bz(ix,iy-1,iz))
     +                     -difz(iz)*(by(ix,iy,iz+1)-by(ix,iy,iz-1))
        feldy(ix,iy,iz) = difz(iz)*(bx(ix,iy,iz+1)-bx(ix,iy,iz-1))
     +                     -difx(ix)*(bz(ix+1,iy,iz)-bz(ix-1,iy,iz))
        feldz(ix,iy,iz) = difx(ix)*(by(ix+1,iy,iz)-by(ix-1,iy,iz))
     +                     -dify(iy)*(bx(ix,iy+1,iz)-bx(ix,iy-1,iz))
        hilfx(ix,iy,iz) = ( sx(ix,iy,iz)*difx(ix)*(
     +                        help(ix+1,iy,iz)-help(ix-1,iy,iz))
     +                 + sy(ix,iy,iz)*dify(iy)*(
     +                        help(ix,iy+1,iz)-help(ix,iy-1,iz))
     +                 + sz(ix,iy,iz)*difz(iz)*(
     +                        help(ix,iy,iz+1)-help(ix,iy,iz-1)) )
     +                           /rho(ix,iy,iz)
   50 continue
c
c    help = .5*pressure;  feld = current density
c
      do 100 iz = 2,nz-1
      do 100 iy = 2,ny-1
      do 100 ix = 2,nx-1
        hilfy(ix,iy,iz) = ( sx(ix,iy,iz)*( feldy(ix,iy,iz)*bz(ix,iy,iz)
     +                              -feldz(ix,iy,iz)*by(ix,iy,iz) )
     +                 + sy(ix,iy,iz)*( feldz(ix,iy,iz)*bx(ix,iy,iz)
     +                              -feldx(ix,iy,iz)*bz(ix,iy,iz) )
     +                 + sz(ix,iy,iz)*( feldx(ix,iy,iz)*by(ix,iy,iz)
     +                              -feldy(ix,iy,iz)*bx(ix,iy,iz) )
     +                            )/rho(ix,iy,iz)
        hilfz(ix,iy,iz) = res(ix,iy,iz)*( 
     +                        feldx(ix,iy,iz)*feldx(ix,iy,iz)
     +                      + feldy(ix,iy,iz)*feldy(ix,iy,iz)
     +                      + feldz(ix,iy,iz)*feldz(ix,iy,iz) )
  100 continue
c
c        hilfx = velocity. * grad (pressure/2)
c        hilfy = currentd. * (vel. x magnetic field)
c        hilfz = res. * currentd.**2
c
      do 200 iz = 2,nz1
      do 200 iy = 2,ny1
      do 200 ix = 2,nx1
	 mass(ldiag)  = mass(ldiag) + hilf(ix,iy,iz)*rho(ix,iy,iz)
	 enmag(ldiag) = enmag(ldiag) + hilf(ix,iy,iz)*heb(ix,iy,iz)
	 enkin(ldiag) = enkin(ldiag) + hilf(ix,iy,iz)*hev(ix,iy,iz)
	 enthe(ldiag) = enthe(ldiag) + hilf(ix,iy,iz)*help(ix,iy,iz)
         vgradp(ldiag)= vgradp(ldiag) + hilf(ix,iy,iz)*hilfx(ix,iy,iz)
         vdjxb(ldiag) = vdjxb(ldiag) + hilf(ix,iy,iz)*hilfy(ix,iy,iz)
         resj2(ldiag) = resj2(ldiag) + hilf(ix,iy,iz)*hilfz(ix,iy,iz)
  200 continue


      ekinpu(ldiag) = -vgradp(ldiag) + vdjxb(ldiag)
      ethpu(ldiag)  =  vgradp(ldiag) + resj2(ldiag)
      ebpu(ldiag)   = -vdjxb(ldiag) - resj2(ldiag)
      if (gamma .ne. 1.0) enthe(ldiag) = enthe(ldiag)/(gamma-1)
c
c    enmag = feldenergy
c    enkin = kinetic energy
c    enthe = thermal energy
c    vgradp = v . grad (p/2)
c    vdjxb =  v . (v x b)
c    resj2 =  res j**2
c    ekinpu  = -vgradp + vdjxb
c    ethpu   =  vgradp + resj2
c    ebpu    = -vdjxb - resj2 = - e mal j
c
      do 600 iz = 2,nz1
      do 600 iy = 2,ny1
      do 600 ix = nxe,nx1
	  massu(ldiag)  = massu(ldiag) + hilf(ix,iy,iz)*rho(ix,iy,iz)
  600 continue

c
c
      if ( mod(ldiag,2) .eq. 1 )  then
 
      do 650 icomp = 1,4
        bmax(icomp,1) = -100.0
        bmin(icomp,1) =  100.0
        vmax(icomp,1) = -100.0
        vmin(icomp,1) =  100.0
        jmax(icomp,1) = -100.0
        jmin(icomp,1) =  100.0
        emax(icomp,1) = -100.0
        emin(icomp,1) =  100.0
  650 continue
      jemin(1,1) = 100.0
      jemin(2,1) = 100.0
      jemax(1,1) = -100.0
      jemax(2,1) = -100.0
      prho(1,1) = -100.0
      prho(2,1) = -100.0
      prho(3,1) = -100.0
      prho(4,1) = -100.0
      
      
      do 700 ix = 2,nx1
        fmax(ix) = bx(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
        fmin(ix) = bx(ix,2,2)
        iymin(ix) = 2
        izmin(ix) = 2
  700 continue
      do 702 iz = 2,nz1
      do 702 iy = 3,ny1
      do 702 ix = 2,nx1
        if (bx(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = bx(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
        if (bx(ix,iy,iz).lt.fmin(ix)) then
         fmin(ix) = bx(ix,iy,iz)
         iymin(ix) = iy
         izmin(ix) = iz
        endif
  702 continue
      do 704 ix = 2,nx1
        if (fmax(ix).gt.bmax(1,1)) then
          bmax(1,1) = fmax(ix)
          bmax(1,2) = x(ix)
          bmax(1,3) = y(iymax(ix))
          bmax(1,4) = z(izmax(ix))
         endif
        if (fmin(ix).lt.bmin(1,1)) then
          bmin(1,1) = fmin(ix)
          bmin(1,2) = x(ix)
          bmin(1,3) = y(iymin(ix))
          bmin(1,4) = z(izmin(ix))
         endif
  704 continue
      do 710 ix = 2,nx1
        fmax(ix) = by(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
        fmin(ix) = by(ix,2,2)
        iymin(ix) = 2
        izmin(ix) = 2
  710 continue
      do 712 iz = 2,nz1
      do 712 iy = 3,ny1
      do 712 ix = 2,nx1
        if (by(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = by(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
        if (by(ix,iy,iz).lt.fmin(ix)) then
         fmin(ix) = by(ix,iy,iz)
         iymin(ix) = iy
         izmin(ix) = iz
        endif
  712 continue
      do 714 ix = 2,nx1
        if (fmax(ix).gt.bmax(2,1)) then
          bmax(2,1) = fmax(ix)
          bmax(2,2) = x(ix)
          bmax(2,3) = y(iymax(ix))
          bmax(2,4) = z(izmax(ix))
         endif
        if (fmin(ix).lt.bmin(2,1)) then
          bmin(2,1) = fmin(ix)
          bmin(2,2) = x(ix)
          bmin(2,3) = y(iymin(ix))
          bmin(2,4) = z(izmin(ix))
         endif
  714 continue
      do 720 ix = 2,nx1
        fmax(ix) = bz(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
        fmin(ix) = bz(ix,2,2)
        iymin(ix) = 2
        izmin(ix) = 2
  720 continue
      do 722 iz = 2,nz1
      do 722 iy = 3,ny1
      do 722 ix = 2,nx1
        if (bz(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = bz(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
        if (bz(ix,iy,iz).lt.fmin(ix)) then
         fmin(ix) = bz(ix,iy,iz)
         iymin(ix) = iy
         izmin(ix) = iz
        endif
  722 continue
      do 724 ix = 2,nx1
        if (fmax(ix).gt.bmax(3,1)) then
          bmax(3,1) = fmax(ix)
          bmax(3,2) = x(ix)
          bmax(3,3) = y(iymax(ix))
          bmax(3,4) = z(izmax(ix))
         endif
        if (fmin(ix).lt.bmin(3,1)) then
          bmin(3,1) = fmin(ix)
          bmin(3,2) = x(ix)
          bmin(3,3) = y(iymin(ix))
          bmin(3,4) = z(izmin(ix))
         endif
  724 continue
      
      do 730 iz = 2,nz1
      do 730 iy = 2,ny1
      do 730 ix = 2,nx1
        hilfx(ix,iy,iz) = sx(ix,iy,iz)/rho(ix,iy,iz)
        hilfy(ix,iy,iz) = sy(ix,iy,iz)/rho(ix,iy,iz)
        hilfz(ix,iy,iz) = sz(ix,iy,iz)/rho(ix,iy,iz)
  730 continue
      do 740 ix = 2,nx1
        fmax(ix) = hilfx(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
        fmin(ix) = hilfx(ix,2,2)
        iymin(ix) = 2
        izmin(ix) = 2
  740 continue
      do 742 iz = 2,nz1
      do 742 iy = 3,ny1
      do 742 ix = 2,nx1
        if (hilfx(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = hilfx(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
        if (hilfx(ix,iy,iz).lt.fmin(ix)) then
         fmin(ix) = hilfx(ix,iy,iz)
         iymin(ix) = iy
         izmin(ix) = iz
        endif
  742 continue
      do 744 ix = 2,nx1
        if (fmax(ix).gt.vmax(1,1)) then
          vmax(1,1) = fmax(ix)
          vmax(1,2) = x(ix)
          vmax(1,3) = y(iymax(ix))
          vmax(1,4) = z(izmax(ix))
         endif
        if (fmin(ix).lt.vmin(1,1)) then
          vmin(1,1) = fmin(ix)
          vmin(1,2) = x(ix)
          vmin(1,3) = y(iymin(ix))
          vmin(1,4) = z(izmin(ix))
         endif
  744 continue
      do 750 ix = 2,nx1
        fmax(ix) = hilfy(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
        fmin(ix) = hilfy(ix,2,2)
        iymin(ix) = 2
        izmin(ix) = 2
  750 continue
      do 752 iz = 2,nz1
      do 752 iy = 3,ny1
      do 752 ix = 2,nx1
        if (hilfy(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = hilfy(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
        if (hilfy(ix,iy,iz).lt.fmin(ix)) then
         fmin(ix) = hilfy(ix,iy,iz)
         iymin(ix) = iy
         izmin(ix) = iz
        endif
  752 continue
      do 754 ix = 2,nx1
        if (fmax(ix).gt.vmax(2,1)) then
          vmax(2,1) = fmax(ix)
          vmax(2,2) = x(ix)
          vmax(2,3) = y(iymax(ix))
          vmax(2,4) = z(izmax(ix))
         endif
        if (fmin(ix).lt.vmin(2,1)) then
          vmin(2,1) = fmin(ix)
          vmin(2,2) = x(ix)
          vmin(2,3) = y(iymin(ix))
          vmin(2,4) = z(izmin(ix))
         endif
  754 continue
      do 760 ix = 2,nx1
        fmax(ix) = hilfz(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
        fmin(ix) = hilfz(ix,2,2)
        iymin(ix) = 2
        izmin(ix) = 2
  760 continue
      do 762 iz = 2,nz1
      do 762 iy = 3,ny1
      do 762 ix = 2,nx1
        if (hilfz(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = hilfz(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
        if (hilfz(ix,iy,iz).lt.fmin(ix)) then
         fmin(ix) = hilfz(ix,iy,iz)
         iymin(ix) = iy
         izmin(ix) = iz
        endif
  762 continue
      do 764 ix = 2,nx1
        if (fmax(ix).gt.vmax(3,1)) then
          vmax(3,1) = fmax(ix)
          vmax(3,2) = x(ix)
          vmax(3,3) = y(iymax(ix))
          vmax(3,4) = z(izmax(ix))
         endif
        if (fmin(ix).lt.vmin(3,1)) then
          vmin(3,1) = fmin(ix)
          vmin(3,2) = x(ix)
          vmin(3,3) = y(iymin(ix))
          vmin(3,4) = z(izmin(ix))
         endif
  764 continue
      
      do 770 iz = 2,nz1
      do 770 iy = 2,ny1
      do 770 ix = 2,nx1
        hilfx(ix,iy,iz) = sqrt( 2.*heb(ix,iy,iz) )
        hilfy(ix,iy,iz) = sqrt( 2.*hev(ix,iy,iz)/rho(ix,iy,iz) )
  770 continue
      do 780 ix = 2,nx1
        fmax(ix) = hilfx(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
  780 continue
      do 782 iz = 2,nz1
      do 782 iy = 3,ny1
      do 782 ix = 2,nx1
        if (hilfx(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = hilfx(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
  782 continue
      do 784 ix = 2,nx1
        if (fmax(ix).gt.bmax(4,1)) then
          bmax(4,1) = fmax(ix)
          bmax(4,2) = x(ix)
          bmax(4,3) = y(iymax(ix))
          bmax(4,4) = z(izmax(ix))
         endif
  784 continue
      do 790 ix = 2,nx1
        fmax(ix) = hilfy(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
  790 continue
      do 792 iz = 2,nz1
      do 792 iy = 3,ny1
      do 792 ix = 2,nx1
        if (hilfy(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = hilfy(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
  792 continue
      do 794 ix = 2,nx1
        if (fmax(ix).gt.vmax(4,1)) then
          vmax(4,1) = fmax(ix)
          vmax(4,2) = x(ix)
          vmax(4,3) = y(iymax(ix))
          vmax(4,4) = z(izmax(ix))
         endif
  794 continue
      
      do 800 iz = 2,nz1
      do 800 iy = 2,ny1
      do 800 ix = 2,nx1
        hilf(ix,iy,iz) = sqrt( feldx(ix,iy,iz)*feldx(ix,iy,iz)
     +                   + feldy(ix,iy,iz)*feldy(ix,iy,iz)
     +                   + feldz(ix,iy,iz)*feldz(ix,iy,iz) )
  800 continue
      do 810 ix = 2,nx1
        fmax(ix) = feldx(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
        fmin(ix) = feldx(ix,2,2)
        iymin(ix) = 2
        izmin(ix) = 2
  810 continue
      do 812 iz = 2,nz1
      do 812 iy = 3,ny1
      do 812 ix = 2,nx1
        if (feldx(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = feldx(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
        if (feldx(ix,iy,iz).lt.fmin(ix)) then
         fmin(ix) = feldx(ix,iy,iz)
         iymin(ix) = iy
         izmin(ix) = iz
        endif
  812 continue
      do 814 ix = 2,nx1
        if (fmax(ix).gt.jmax(1,1)) then
          jmax(1,1) = fmax(ix)
          jmax(1,2) = x(ix)
          jmax(1,3) = y(iymax(ix))
          jmax(1,4) = z(izmax(ix))
         endif
        if (fmin(ix).lt.jmin(1,1)) then
          jmin(1,1) = fmin(ix)
          jmin(1,2) = x(ix)
          jmin(1,3) = y(iymin(ix))
          jmin(1,4) = z(izmin(ix))
         endif
  814 continue
      do 820 ix = 2,nx1
        fmax(ix) = feldy(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
        fmin(ix) = feldy(ix,2,2)
        iymin(ix) = 2
        izmin(ix) = 2
  820 continue
      do 822 iz = 2,nz1
      do 822 iy = 3,ny1
      do 822 ix = 2,nx1
        if (feldy(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = feldy(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
        if (feldy(ix,iy,iz).lt.fmin(ix)) then
         fmin(ix) = feldy(ix,iy,iz)
         iymin(ix) = iy
         izmin(ix) = iz
        endif
  822 continue
      do 824 ix = 2,nx1
        if (fmax(ix).gt.jmax(2,1)) then
          jmax(2,1) = fmax(ix)
          jmax(2,2) = x(ix)
          jmax(2,3) = y(iymax(ix))
          jmax(2,4) = z(izmax(ix))
         endif
        if (fmin(ix).lt.jmin(2,1)) then
          jmin(2,1) = fmin(ix)
          jmin(2,2) = x(ix)
          jmin(2,3) = y(iymin(ix))
          jmin(2,4) = z(izmin(ix))
         endif
  824 continue
      do 830 ix = 2,nx1
        fmax(ix) = feldz(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
        fmin(ix) = feldz(ix,2,2)
        iymin(ix) = 2
        izmin(ix) = 2
  830 continue
      do 832 iz = 2,nz1
      do 832 iy = 3,ny1
      do 832 ix = 2,nx1
        if (feldz(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = feldz(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
        if (feldz(ix,iy,iz).lt.fmin(ix)) then
         fmin(ix) = feldz(ix,iy,iz)
         iymin(ix) = iy
         izmin(ix) = iz
        endif
  832 continue
      do 834 ix = 2,nx1
        if (fmax(ix).gt.jmax(3,1)) then
          jmax(3,1) = fmax(ix)
          jmax(3,2) = x(ix)
          jmax(3,3) = y(iymax(ix))
          jmax(3,4) = z(izmax(ix))
         endif
        if (fmin(ix).lt.jmin(3,1)) then
          jmin(3,1) = fmin(ix)
          jmin(3,2) = x(ix)
          jmin(3,3) = y(iymin(ix))
          jmin(3,4) = z(izmin(ix))
         endif
  834 continue
      do 840 ix = 2,nx1
        fmax(ix) = hilf(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
  840 continue
      do 842 iz = 2,nz1
      do 842 iy = 3,ny1
      do 842 ix = 2,nx1
        if (hilf(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = hilf(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
  842 continue
      do 844 ix = 2,nx1
        if (fmax(ix).gt.jmax(4,1)) then
          jmax(4,1) = fmax(ix)
          jmax(4,2) = x(ix)
          jmax(4,3) = y(iymax(ix))
          jmax(4,4) = z(izmax(ix))
         endif
  844 continue

      do 850 iz = 2,nz1
      do 850 iy = 2,ny1
      do 850 ix = 2,nx1
        hilfx(ix,iy,iz) = res(ix,iy,iz)*feldx(ix,iy,iz) - (
     +                     sy(ix,iy,iz)*bz(ix,iy,iz)
     +                   - sz(ix,iy,iz)*by(ix,iy,iz) )/rho(ix,iy,iz)
        hilfy(ix,iy,iz) = res(ix,iy,iz)*feldy(ix,iy,iz) - (
     +                     sz(ix,iy,iz)*bx(ix,iy,iz)
     +                   - sx(ix,iy,iz)*bz(ix,iy,iz) )/rho(ix,iy,iz)
        hilfz(ix,iy,iz) = res(ix,iy,iz)*feldz(ix,iy,iz) - (
     +                     sx(ix,iy,iz)*by(ix,iy,iz)
     +                   - sy(ix,iy,iz)*bx(ix,iy,iz) )/rho(ix,iy,iz)
  850 continue
      do 855 iz = 2,nz1
      do 855 iy = 2,ny1
      do 855 ix = 2,nx1
        hilf(ix,iy,iz) = sqrt( hilfx(ix,iy,iz)*hilfx(ix,iy,iz)
     +                   + hilfy(ix,iy,iz)*hilfy(ix,iy,iz)
     +                   + hilfz(ix,iy,iz)*hilfz(ix,iy,iz) )
  855 continue
      do 860 ix = 2,nx1
        fmax(ix) = hilfx(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
        fmin(ix) = hilfx(ix,2,2)
        iymin(ix) = 2
        izmin(ix) = 2
  860 continue
      do 862 iz = 2,nz1
      do 862 iy = 3,ny1
      do 862 ix = 2,nx1
        if (hilfx(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = hilfx(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
        if (hilfx(ix,iy,iz).lt.fmin(ix)) then
         fmin(ix) = hilfx(ix,iy,iz)
         iymin(ix) = iy
         izmin(ix) = iz
        endif
  862 continue
      do 864 ix = 2,nx1
        if (fmax(ix).gt.emax(1,1)) then
          emax(1,1) = fmax(ix)
          emax(1,2) = x(ix)
          emax(1,3) = y(iymax(ix))
          emax(1,4) = z(izmax(ix))
         endif
        if (fmin(ix).lt.emin(1,1)) then
          emin(1,1) = fmin(ix)
          emin(1,2) = x(ix)
          emin(1,3) = y(iymin(ix))
          emin(1,4) = z(izmin(ix))
         endif
  864 continue
      do 870 ix = 2,nx1
        fmax(ix) = hilfy(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
        fmin(ix) = hilfy(ix,2,2)
        iymin(ix) = 2
        izmin(ix) = 2
  870 continue
      do 872 iz = 2,nz1
      do 872 iy = 3,ny1
      do 872 ix = 2,nx1
        if (hilfy(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = hilfy(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
        if (hilfy(ix,iy,iz).lt.fmin(ix)) then
         fmin(ix) = hilfy(ix,iy,iz)
         iymin(ix) = iy
         izmin(ix) = iz
        endif
  872 continue
      do 874 ix = 2,nx1
        if (fmax(ix).gt.emax(2,1)) then
          emax(2,1) = fmax(ix)
          emax(2,2) = x(ix)
          emax(2,3) = y(iymax(ix))
          emax(2,4) = z(izmax(ix))
         endif
        if (fmin(ix).lt.emin(2,1)) then
          emin(2,1) = fmin(ix)
          emin(2,2) = x(ix)
          emin(2,3) = y(iymin(ix))
          emin(2,4) = z(izmin(ix))
         endif
  874 continue
      do 880 ix = 2,nx1
        fmax(ix) = hilfz(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
        fmin(ix) = hilfz(ix,2,2)
        iymin(ix) = 2
        izmin(ix) = 2
  880 continue
      do 882 iz = 2,nz1
      do 882 iy = 3,ny1
      do 882 ix = 2,nx1
        if (hilfz(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = hilfz(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
        if (hilfz(ix,iy,iz).lt.fmin(ix)) then
         fmin(ix) = hilfz(ix,iy,iz)
         iymin(ix) = iy
         izmin(ix) = iz
        endif
  882 continue
      do 884 ix = 2,nx1
        if (fmax(ix).gt.emax(3,1)) then
          emax(3,1) = fmax(ix)
          emax(3,2) = x(ix)
          emax(3,3) = y(iymax(ix))
          emax(3,4) = z(izmax(ix))
         endif
        if (fmin(ix).lt.emin(3,1)) then
          emin(3,1) = fmin(ix)
          emin(3,2) = x(ix)
          emin(3,3) = y(iymin(ix))
          emin(3,4) = z(izmin(ix))
         endif
  884 continue
      do 890 ix = 2,nx1
        fmax(ix) = hilf(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
  890 continue
      do 892 iz = 2,nz1
      do 892 iy = 3,ny1
      do 892 ix = 2,nx1
        if (hilf(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = hilf(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
  892 continue
      do 894 ix = 2,nx1
        if (fmax(ix).gt.emax(4,1)) then
          emax(4,1) = fmax(ix)
          emax(4,2) = x(ix)
          emax(4,3) = y(iymax(ix))
          emax(4,4) = z(izmax(ix))
         endif
  894 continue
  
      do 900 iz = 2,nz1
      do 900 iy = 2,ny1
      do 900 ix = 2,nx1
        if (heb(ix,iy,iz).ne.0.0) then
           hilf(ix,iy,iz) = ( bx(ix,iy,iz)*feldx(ix,iy,iz)
     +                       + by(ix,iy,iz)*feldy(ix,iy,iz)
     +                       + bz(ix,iy,iz)*feldz(ix,iy,iz) )
     +                       /sqrt(2.*heb(ix,iy,iz))
        else
            hilf(ix,iy,iz) = 0.0
        end if  
  900 continue
      do 905 iz = 2,nz1
      do 905 iy = 2,ny1
      do 905 ix = 2,nx1
        hilfx(ix,iy,iz)=hilf(ix,iy,iz)*res(ix,iy,iz)
  905 continue
      do 910 ix = 2,nx1
        fmax(ix) = hilf(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
        fmin(ix) = hilf(ix,2,2)
        iymin(ix) = 2
        izmin(ix) = 2
  910 continue
      do 912 iz = 2,nz1
      do 912 iy = 3,ny1
      do 912 ix = 2,nx1
        if (hilf(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = hilf(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
        if (hilf(ix,iy,iz).lt.fmin(ix)) then
         fmin(ix) = hilf(ix,iy,iz)
         iymin(ix) = iy
         izmin(ix) = iz
        endif
  912 continue
      do 914 ix = 2,nx1
        if (fmax(ix).gt.jemax(1,1)) then
          jemax(1,1) = fmax(ix)
          jemax(1,2) = x(ix)
          jemax(1,3) = y(iymax(ix))
          jemax(1,4) = z(izmax(ix))
         endif
        if (fmin(ix).lt.jemin(1,1)) then
          jemin(1,1) = fmin(ix)
          jemin(1,2) = x(ix)
          jemin(1,3) = y(iymin(ix))
          jemin(1,4) = z(izmin(ix))
         endif
  914 continue
      do 920 ix = 2,nx1
        fmax(ix) = hilfx(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
        fmin(ix) = hilfx(ix,2,2)
        iymin(ix) = 2
        izmin(ix) = 2
  920 continue
      do 922 iz = 2,nz1
      do 922 iy = 3,ny1
      do 922 ix = 2,nx1
        if (hilfx(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = hilfx(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
        if (hilfx(ix,iy,iz).lt.fmin(ix)) then
         fmin(ix) = hilfx(ix,iy,iz)
         iymin(ix) = iy
         izmin(ix) = iz
        endif
  922 continue
      do 924 ix = 2,nx1
        if (fmax(ix).gt.jemax(2,1)) then
          jemax(2,1) = fmax(ix)
          jemax(2,2) = x(ix)
          jemax(2,3) = y(iymax(ix))
          jemax(2,4) = z(izmax(ix))
         endif
        if (fmin(ix).lt.jemin(2,1)) then
          jemin(2,1) = fmin(ix)
          jemin(2,2) = x(ix)
          jemin(2,3) = y(iymin(ix))
          jemin(2,4) = z(izmin(ix))
         endif
  924 continue
  
      do 930 iz = 2,nz1
      do 930 iy = 2,ny1
      do 930 ix = 2,nx1
        hilfx(ix,iy,iz)=2.*help(ix,iy,iz)+2.*heb(ix,iy,iz)
        hilfy(ix,iy,iz)=2.*help(ix,iy,iz)
  930 continue
      do 940 ix = 2,nx1
        fmax(ix) = hilfx(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
  940 continue
      do 942 iz = 2,nz1
      do 942 iy = 3,ny1
      do 942 ix = 2,nx1
        if (hilfx(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = hilfx(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
  942 continue
      do 944 ix = 2,nx1
        if (fmax(ix).gt.prho(1,1)) then
          prho(1,1) = fmax(ix)
          prho(1,2) = x(ix)
          prho(1,3) = y(iymax(ix))
          prho(1,4) = z(izmax(ix))
         endif
  944 continue
      do 950 ix = 2,nx1
        fmax(ix) = hilfy(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
  950 continue
      do 952 iz = 2,nz1
      do 952 iy = 3,ny1
      do 952 ix = 2,nx1
        if (hilfy(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = hilfy(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
  952 continue
      do 954 ix = 2,nx1
        if (fmax(ix).gt.prho(2,1)) then
          prho(2,1) = fmax(ix)
          prho(2,2) = x(ix)
          prho(2,3) = y(iymax(ix))
          prho(2,4) = z(izmax(ix))
         endif
  954 continue
      do 960 ix = 2,nx1
        fmax(ix) = rho(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
  960 continue
      do 962 iz = 2,nz1
      do 962 iy = 3,ny1
      do 962 ix = 2,nx1
        if (rho(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = rho(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
  962 continue
      do 964 ix = 2,nx1
        if (fmax(ix).gt.prho(3,1)) then
          prho(3,1) = fmax(ix)
          prho(3,2) = x(ix)
          prho(3,3) = y(iymax(ix))
          prho(3,4) = z(izmax(ix))
         endif
  964 continue
      do 970 ix = 2,nx1
        fmax(ix) = res(ix,2,2)
        iymax(ix) = 2
        izmax(ix) = 2
  970 continue
      do 972 iz = 2,nz1
      do 972 iy = 3,ny1
      do 972 ix = 2,nx1
        if (res(ix,iy,iz).gt.fmax(ix)) then
         fmax(ix) = res(ix,iy,iz)
         iymax(ix) = iy
         izmax(ix) = iz
        endif
  972 continue
      do 974 ix = 2,nx1
        if (fmax(ix).gt.prho(4,1)) then
          prho(4,1) = fmax(ix)
          prho(4,2) = x(ix)
          prho(4,3) = y(iymax(ix))
          prho(4,4) = z(izmax(ix))
         endif
  974 continue

      write(25,11) time
      write(25,12) bmax(1,1),bmax(1,2),bmax(1,3),bmax(1,4),
     +             bmin(1,1),bmin(1,2),bmin(1,3),bmin(1,4)
      write(25,13) bmax(2,1),bmax(2,2),bmax(2,3),bmax(2,4),
     +             bmin(2,1),bmin(2,2),bmin(2,3),bmin(2,4)
      write(25,13) bmax(3,1),bmax(3,2),bmax(3,3),bmax(3,4),
     +             bmin(3,1),bmin(3,2),bmin(3,3),bmin(3,4)
      write(25,13) bmax(4,1),bmax(4,2),bmax(4,3),bmax(4,4)
      write(25,14) vmax(1,1),vmax(1,2),vmax(1,3),vmax(1,4),
     +             vmin(1,1),vmin(1,2),vmin(1,3),vmin(1,4)
      write(25,13) vmax(2,1),vmax(2,2),vmax(2,3),vmax(2,4),
     +             vmin(2,1),vmin(2,2),vmin(2,3),vmin(2,4)
      write(25,13) vmax(3,1),vmax(3,2),vmax(3,3),vmax(3,4),
     +             vmin(3,1),vmin(3,2),vmin(3,3),vmin(3,4)
      write(25,13) vmax(4,1),vmax(4,2),vmax(4,3),vmax(4,4)
      write(25,15) jmax(1,1),jmax(1,2),jmax(1,3),jmax(1,4),
     +             jmin(1,1),jmin(1,2),jmin(1,3),jmin(1,4)
      write(25,13) jmax(2,1),jmax(2,2),jmax(2,3),jmax(2,4),
     +             jmin(2,1),jmin(2,2),jmin(2,3),jmin(2,4)
      write(25,13) jmax(3,1),jmax(3,2),jmax(3,3),jmax(3,4),
     +             jmin(3,1),jmin(3,2),jmin(3,3),jmin(3,4)
      write(25,13) jmax(4,1),jmax(4,2),jmax(4,3),jmax(4,4)
      write(25,16) emax(1,1),emax(1,2),emax(1,3),emax(1,4),
     +             emin(1,1),emin(1,2),emin(1,3),emin(1,4)
      write(25,13) emax(2,1),emax(2,2),emax(2,3),emax(2,4),
     +             emin(2,1),emin(2,2),emin(2,3),emin(2,4)
      write(25,13) emax(3,1),emax(3,2),emax(3,3),emax(3,4),
     +             emin(3,1),emin(3,2),emin(3,3),emin(3,4)
      write(25,13) emax(4,1),emax(4,2),emax(4,3),emax(4,4)
      write(25,17) jemax(1,1),jemax(1,2),jemax(1,3),jemax(1,4),
     +             jemin(1,1),jemin(1,2),jemin(1,3),jemin(1,4)
      write(25,13) jemax(2,1),jemax(2,2),jemax(2,3),jemax(2,4),
     +             jemin(2,1),jemin(2,2),jemin(2,3),jemin(2,4)
      write(25,18) prho(1,1),prho(1,2),prho(1,3),prho(1,4)
      write(25,13) prho(2,1),prho(2,2),prho(2,3),prho(2,4)
      write(25,13) prho(3,1),prho(3,2),prho(3,3),prho(3,4)
      write(25,19) prho(4,1),prho(4,2),prho(4,3),prho(4,4)

   11 format(f8.2)
   12 format(f11.4,3f8.2,f11.4,3f8.2,3x,'B')
   13 format(f11.4,3f8.2,f11.4,3f8.2)
   14 format(f11.4,3f8.2,f11.4,3f8.2,3x,'V')
   15 format(f11.4,3f8.2,f11.4,3f8.2,3x,'J')
   16 format(f11.4,3f8.2,f11.4,3f8.2,3x,'E')
   17 format(f11.4,3f8.2,f11.4,3f8.2,3x,'J,E ll')
   18 format(f11.4,3f8.2,36x,'PB,P,Rho')
   19 format(f11.4,3f8.2,36x,'Res')
      end if
c
      return
      end
	subroutine diag2
c***********************************************************************
	include 'misflin'
c
        integer  iz,iy,nxe
        real     hx2(ny,nz),hxy1(nz),hxz1(ny)
c.....................................................
c   flux of mass, energies and magnetic field through
c   x = const (here x=0)
c.....................................................
      if (zentr(1)) then
        nxe = (nx+1)/2
      else
        nxe = 2
      end if

      do 20 iz = 1,nz
      do 20 iy = 1,ny
	hx2(iy,iz)  = 0.25/dify(iy)/difz(iz)
	if (iy .eq. 2 .or. iy .eq. ny1) hx2(iy,iz) = 0.5*hx2(iy,iz) 
	if (iz .eq. 2 .or. iz .eq. nz1) hx2(iy,iz) = 0.5*hx2(iy,iz) 
   20 continue
      do 40 iz = 1,nz
      do 40 iy = 1,ny
        help(nxe,iy,iz)  = u(nxe,iy,iz)**gamma
	hilfx(nxe,iy,iz) = sx(nxe,iy,iz)/rho(nxe,iy,iz)
   40 continue

      do 50 iz = 2,nz1
      do 50 iy = 2,ny1
        feldx(nxe,iy,iz) =  hilfx(nxe,iy,iz)*( 2.*heb(nxe,iy,iz)
     +                  - bx(nxe,iy,iz)/rho(nxe,iy,iz)*( 
     +                       sx(nxe,iy,iz)*bx(nxe,iy,iz)
     +                     + sy(nxe,iy,iz)*by(nxe,iy,iz)
     +                     + sz(nxe,iy,iz)*bz(nxe,iy,iz) )  )
        feldy(nxe,iy,iz) = difz(iz)*(bx(nxe,iy,iz+1)-bx(nxe,iy,iz-1))
     +                     -difx(nxe)*(bz(nxe+1,iy,iz)-bz(nxe-1,iy,iz))
        feldz(nxe,iy,iz) = difx(nxe)*(by(nxe+1,iy,iz)-by(nxe-1,iy,iz))
     +                     -dify(iy)*(bx(nxe,iy+1,iz)-bx(nxe,iy-1,iz))
   50 continue
c
      do 100 iz = 2,nz1
      do 100 iy = 2,ny1
       hilfz(nxe,iy,iz) = res(nxe,iy,iz)
     +                        *( feldy(nxe,iy,iz)*bz(nxe,iy,iz)
     +                           -feldz(nxe,iy,iz)*by(nxe,iy,iz) )
  100 continue
c
c        hilf  = dely*delz
c        help  = druck/2
c        hilfx = vx
c        hev = 0.5*rho*v**2
c        hilfz = res*(strom x b)_xcomp 
c        feldx = -(v x b) x b 
c        feldx+hilfz = (e x b)_xcomp 
c        feldy = stromy 
c        feldz = stromz 
c
      do 200 iz = 2,nz1
      do 200 iy = 2,ny1
	  frho(ldiag)  = frho(ldiag) - hx2(iy,iz)*sx(nxe,iy,iz)
	  fethe(ldiag) = fethe(ldiag) 
     +                   - hx2(iy,iz)*hilfx(nxe,iy,iz)*help(nxe,iy,iz)
	  fekin(ldiag) = fekin(ldiag)  
     +                   - hx2(iy,iz)*hilfx(nxe,iy,iz)*hev(nxe,iy,iz)
	  febr(ldiag) = febr(ldiag) 
     +                   - hx2(iy,iz)*feldx(nxe,iy,iz)
	  febi(ldiag) = febi(ldiag) 
     +                   - hx2(iy,iz)*hilfz(nxe,iy,iz)
	  fphi(ldiag)  = fphi(ldiag) - hx2(iy,iz)*bx(nxe,iy,iz)
  200 continue

      do 400 iz = 1,nz
	hxy1(iz)  = 0.5/difz(iz)
	if (iz .eq. 2 .or. iz .eq. nz1) hxy1(iz) = 0.5*hxy1(iz) 
  400 continue

      do 410 iy = 1,ny
	hxz1(iy)  = 0.5/dify(iy)
	if (iy .eq. 2 .or. iy .eq. ny1) hxz1(iy) = 0.5*hxz1(iy) 
  410 continue

      do 420 iy = 2,ny1
        fphipi(ldiag) = fphipi(ldiag) + hxz1(iy)*(
     +             ( sz(nxe,iy,nz1)*bx(nxe,iy,nz1)
     +              - sx(nxe,iy,nz1)*bz(nxe,iy,nz1) )/rho(nxe,iy,nz1)
     +             - ( sz(nxe,iy,2)*bx(nxe,iy,2) 
     +                - sx(nxe,iy,2)*bz(nxe,iy,2) )/rho(nxe,iy,2) )
        fphipr(ldiag) = fphipr(ldiag) + hxz1(iy)*(
     +              - res(nxe,iy,nz1)*feldy(nxe,iy,nz1)
     +                  + res(nxe,iy,2)*feldy(nxe,iy,2) )
  420 continue

      do 430 iz = 2, nz1
        fphipi(ldiag) = fphipi(ldiag) + hxy1(iz)*(
     +            ( sx(nxe,2,iz)*by(nxe,2,iz) 
     +             - sy(nxe,2,iz)*bx(nxe,2,iz) )/rho(nxe,2,iz)
     +         - ( sx(nxe,ny1,iz)*by(nxe,ny1,iz) 
     +            - sy(nxe,ny1,iz)*bx(nxe,ny1,iz) )/rho(nxe,ny1,iz) )
        fphipr(ldiag) = fphipr(ldiag)  + hxy1(iz)*(
     +              - res(nxe,2,iz)*feldz(nxe,2,iz)
     +                  + res(nxe,ny1,iz)*feldz(nxe,ny1,iz) )
  430 continue

      feb(ldiag) = febi(ldiag) + febr(ldiag) 
      fphip(ldiag) = fphipi(ldiag) + fphipr(ldiag)
      if (gamma .ne. 1.0) fethe(ldiag) = fethe(ldiag)*gamma/(gamma-1)

c
c alle groesse bez. fluss in negative x richtung durch x = 0
c
c   flux of 
c    frho = mass
c    fekin = kinetic energie
c    fethe = thermal energie
c    feb = fieldenergie
c      index: i-ideale contrib, r-resistive contrib
c    fphi = magnetic fieldes
c    fphip = flux change of magnetic field
c      index: i-ideale contrib, r-resistive contrib
c
      return
      end
	subroutine diag3
c***********************************************************************
	include 'misflin'
c
        integer  ix,iy,iz,ix0,iz0,stpy,
     +           ixp(ny),izp(ny),
     +           ixbp(ny),izbp(ny),ixbm(ny),izbm(ny)
	real     pmax(ny),bnp(ny),bnm(ny),vx3,vy3,vz3,
     +           xp(ny),zp(ny),xbp(ny),zbp(ny),xbm(ny),zbm(ny),
     +           dxp,dxm,dfpx,dfmx,dzp,dzm,dfpz,dfmz,
     +           hnenx,hnenz,alphx,alphz,betax,betaz,
     +           delx,delz,delb
c.....................................................
c   properties of the reconnection bulge
c.....................................................
      stpy = 4
      do 10 iy = 1,ny
        pmax(iy) = 1.0 + pmsp + 0.01
        xp(iy) = 0.0
        zp(iy) = 0.0
        ixp(iy) = 1
        izp(iy) = 1
        bnp(iy) = 0.001
        xbp(iy) = 0.0
        zbp(iy) = 0.0
        ixbp(iy) = 1
        izbp(iy) = 1
        bnm(iy) = -0.001
        xbm(iy) = 0.0
        zbm(iy) = 0.0
        ixbm(iy) = 1
        izbm(iy) = 1
   10 continue

      do 40 iz = 2, nz1
      do 40 iy = 2, ny1
      do 40 ix = 2, nx1
        help(ix,iy,iz)  = 2.0*(u(ix,iy,iz)**gamma + heb(ix,iy,iz))
   40 continue

      do 100 iy = 2,ny1,stpy
      do 100 iz = 3,nz-2
      do 100 ix = 3,nx-2
        if (help(ix,iy,iz) .gt. pmax(iy)) then
          pmax(iy) = help(ix,iy,iz)
          ixp(iy) = ix
          izp(iy) = iz
        end if
        if (bx(ix,iy,iz) .gt. bnp(iy)) then
          bnp(iy) = bx(ix,iy,iz)
          ixbp(iy) = ix
          izbp(iy) = iz
        end if
        if (bx(ix,iy,iz) .lt. bnm(iy)) then
          bnm(iy) = bx(ix,iy,iz)
          ixbm(iy) = ix
          izbm(iy) = iz
        end if
  100 continue

      do 200 iy = 2,ny1,stpy
       if ( izp(iy).lt. nz .and. izp(iy) .gt. 1 .and.
     +      ixp(iy).lt. nx1 .and. ixp(iy) .gt. 2 ) then
        ix0 = ixp(iy)
        iz0 = izp(iy)
        dxp = x(ix0+1)-x(ix0)
        dxm = x(ix0)-x(ix0-1)
        dfpx = help(ix0+1,iy,iz0)-help(ix0,iy,iz0)
        dfmx = help(ix0-1,iy,iz0)-help(ix0,iy,iz0)
        dzp = z(iz0+1)-z(iz0)
        dzm = z(iz0)-z(iz0-1)
        dfpz = help(ix0,iy,iz0+1)-help(ix0,iy,iz0)
        dfmz = help(ix0,iy,iz0-1)-help(ix0,iy,iz0)
        hnenx = dxp*dxm*(dxp+dxm)
        hnenz = dzp*dzm*(dzp+dzm)
        alphx = (dfpx*dxm+dfmx*dxp)/hnenx
        alphz = (dfpz*dzm+dfmz*dzp)/hnenz
        betax = (dfpx*dxm*dxm-dfmx*dxp*dxp)/hnenx
        betaz = (dfpz*dzm*dzm-dfmz*dzp*dzp)/hnenz
        if (alphx.ne.0) then
          delx = -0.5*betax/alphx
        else 
          delx = 0.0
        end if
        if (alphz.ne.0) then
          delz = -0.5*betaz/alphz
        else 
          delz = 0.0
        end if
        xp(iy) = x(ix0) + delx
        zp(iy) = z(iz0) + delz
        pmax(iy) = pmax(iy) + alphx*delx*delx + betax*delx
     +                      + alphz*delz*delz + betaz*delz
       else
        xp(iy) = 0
        zp(iy) = 0
        pmax(iy) = 0
      end if
  200 continue

      do 300 iy = 2,ny1,stpy
       if ( izbp(iy).lt. nz .and. izbp(iy) .gt. 1 .and.
     +      ixbp(iy).lt. nx1 .and. ixbp(iy) .gt. 2 ) then
        ix0 = ixbp(iy)
        iz0 = izbp(iy)
        dxp = x(ix0+1)-x(ix0)
        dxm = x(ix0)-x(ix0-1)
        dfpx = bx(ix0+1,iy,iz0)-bx(ix0,iy,iz0)
        dfmx = bx(ix0-1,iy,iz0)-bx(ix0,iy,iz0)
        dzp = z(iz0+1)-z(iz0)
        dzm = z(iz0)-z(iz0-1)
        dfpz = bx(ix0,iy,iz0+1)-bx(ix0,iy,iz0)
        dfmz = bx(ix0,iy,iz0-1)-bx(ix0,iy,iz0)
        hnenx = dxp*dxm*(dxp+dxm)
        hnenz = dzp*dzm*(dzp+dzm)
        alphx = (dfpx*dxm+dfmx*dxp)/hnenx
        alphz = (dfpz*dzm+dfmz*dzp)/hnenz
        betax = (dfpx*dxm*dxm-dfmx*dxp*dxp)/hnenx
        betaz = (dfpz*dzm*dzm-dfmz*dzp*dzp)/hnenz
        if (alphx.ne.0) then
          delx = -0.5*betax/alphx
        else 
          delx = 0.0
        end if
        if (alphz.ne.0) then
          delz = -0.5*betaz/alphz
        else 
          delz = 0.0
        end if
        xbp(iy) = x(ix0) + delx
        zbp(iy) = z(iz0) + delz
        bnp(iy) = bnp(iy) + alphx*delx*delx + betax*delx
     +                      + alphz*delz*delz + betaz*delz
       else
        xbp(iy) = 0
        zbp(iy) = 0
        bnp(iy) = 0
      end if
  300 continue

      do 400 iy = 2,ny1,stpy
       if ( izbm(iy).lt. nz .and. izbm(iy) .gt. 1 .and.
     +      ixbm(iy).lt. nx1 .and. ixbm(iy) .gt. 2 ) then
        ix0 = ixbm(iy)
        iz0 = izbm(iy)
        dxp = x(ix0+1)-x(ix0)
        dxm = x(ix0)-x(ix0-1)
        dfpx = bx(ix0+1,iy,iz0)-bx(ix0,iy,iz0)
        dfmx = bx(ix0-1,iy,iz0)-bx(ix0,iy,iz0)
        dzp = z(iz0+1)-z(iz0)
        dzm = z(iz0)-z(iz0-1)
        dfpz = bx(ix0,iy,iz0+1)-bx(ix0,iy,iz0)
        dfmz = bx(ix0,iy,iz0-1)-bx(ix0,iy,iz0)
        hnenx = dxp*dxm*(dxp+dxm)
        hnenz = dzp*dzm*(dzp+dzm)
        alphx = (dfpx*dxm+dfmx*dxp)/hnenx
        alphz = (dfpz*dzm+dfmz*dzp)/hnenz
        betax = (dfpx*dxm*dxm-dfmx*dxp*dxp)/hnenx
        betaz = (dfpz*dzm*dzm-dfmz*dzp*dzp)/hnenz
        if (alphx.ne.0) then
          delx = -0.5*betax/alphx
        else 
          delx = 0.0
        end if
        if (alphz.ne.0) then
          delz = -0.5*betaz/alphz
        else 
          delz = 0.0
        end if
        xbm(iy) = x(ix0) + delx
        zbm(iy) = z(iz0) + delz
        bnm(iy) = bnm(iy) + alphx*delx*delx + betax*delx
     +                      + alphz*delz*delz + betaz*delz
       else
        xbm(iy) = 0
        zbm(iy) = 0
        bnm(iy) = 0
      end if
  400 continue

      write(19,1) time
      write(19,2)
      do 800 iy = 2,ny1,stpy
        ix0=ixp(iy)
        iz0 = izp(iy)
        vx3 = sx(ix0,iy,iz0)/rho(ix0,iy,iz0)
        vy3 = sy(ix0,iy,iz0)/rho(ix0,iy,iz0)
        vz3 = sz(ix0,iy,iz0)/rho(ix0,iy,iz0)
        delb = bnp(iy)-bnm(iy)
        delx = xbp(iy)-xbm(iy)
        delz = zbp(iy)-zbm(iy)
        write(19,3) y(iy),pmax(iy),xp(iy),zp(iy),rho(ix0,iy,iz0),
     +              bx(ix0,iy,iz0),by(ix0,iy,iz0),bz(ix0,iy,iz0)
        write(19,4) bnp(iy),xbp(iy),zbp(iy),
     +              vx3,vy3,vz3
        write(19,5) bnm(iy),xbm(iy),zbm(iy),
     +              delb,delx,delz
  800 continue

    1 format(/1x,'time = ',f8.2)
    2 format(1x,'   y   !  pmax  !   xp   !   zp   !',
     +            '  rho  !   bx   !   by   !   bz   !',
     +      /1x,'--------!  bxmax !  xbmax !  zbmax !',
     +            '------!   vx   !   vy   !   vz   !',
     +      /1x,'---------!  bxmin !  xbmin !  zbmin !',
     +            '!!  difb  !  difx  !  difz  !',
     +      /1x,'-----------------------------------------',
     +                    '------------------------------')
    3 format(1x,f7.3,3f9.3,f8.3,3f9.3)
    4 format(9x,3f9.3,7x,3f9.3)
    5 format(10x,3f9.3,2x,3f9.3)
c
      return
      end
	subroutine diagout
c***********************************************************************
	include 'misflin'
c
        real rrho,rekin,rethe,reb,rebi,rebr
	integer  i
c ....................................................
c     output magdia1(17) and magdia2(18)
c ....................................................
      do 1100 i = 1,ldiag
            rrho = 0.0
            rekin = 0.0
            rethe = 0.0
            reb = 0.0
            rebi = 0.0
            rebr = 0.0
          if (abs(fphi(i)).gt. 0.01) then 
            rrho = frho(i)/fphi(i)
            rekin = fekin(i)/fphi(i)
            rethe = fethe(i)/fphi(i)
            reb = feb(i)/fphi(i)
            rebi = febi(i)/fphi(i)
            rebr = febr(i)/fphi(i)
          end if

	  write(17,1) tdiag(i),mass(i),massu(i),
     +                         enkin(i),enmag(i),enthe(i)
	  write(18,2) tdiag(i),vgradp(i),vdjxb(i),resj2(i),
     +                ekinpu(i),ebpu(i),ethpu(i)
	  write(27,3) tdiag(i),fphi(i),fphip(i),fphipi(i),fphipr(i)
	  write(28,2) tdiag(i),frho(i),fekin(i),fethe(i),
     +                feb(i),febi(i),febr(i)
	  write(28,4) rrho,rekin,rethe,
     +                reb,rebi,rebr
 1100 continue
      do 1200 i = 1,ndiag
	  mass(i)  = 0.
	  massu(i)  = 0.
	  enmag(i) = 0.
 	  enkin(i) = 0.
	  enthe(i) = 0.
	  vgradp(i) = 0.
	  vdjxb(i) = 0.
	  resj2(i) = 0.
	  ekinpu(i)  = 0.
	  ethpu(i)   = 0.
	  ebpu(i)    = 0.
	  tdiag(i)   = 0.
          fphi(i)   = 0.
          fphip(i)  = 0.
          fphipi(i) = 0.
          fphipr(i) = 0.
          frho(i)   = 0.
          fekin(i)  = 0.
          fethe(i)  = 0.
          feb(i)    = 0.
          febi(i)   = 0.
          febr(i)   = 0.
 1200 continue
c
    1 format(1x,f7.3,' ',f12.2,' ',f12.2,' ',f12.3,' ',f12.2,' ',f12.2)
    2 format(1x,f7.3,' ',f10.3,' ',f10.3,' ',f10.3,' ',f10.3,
     +               ' ',f10.3,' ',f10.3 )
    3 format(1x,f7.3,' ',f12.3,' ',f12.3,' ',f12.3,' ',f12.3)
    4 format(5x,    '  ',f10.4,' ',f10.4,' ',f10.4,' ',f10.4,
     +               ' ',f10.4,' ',f10.4 )
c
      return
      end
	subroutine resist
c***********************************************************************
	include 'misflin'
c
	integer  ix,iy,iz,ixanf
	real     eta0,eta1,jcrit1,jcrit2,epar1,ppp1(nx),ppp2(nx)
c
c************************************************
c    3.a resistivity
c    -------------------------------
c************************************************
c
      if (etasw .eq. 1) then
        eta0 = eta*( 1.0-exp(-time/5.0) )
        eta1 = 0.02*eta0
        do 10 iz = 1,nzz
        do 10 iy = 1,ny
        do 10 ix = 1,nx
           res(ix,iy,iz) = eta0*prof(ix,iy,iz) + eta1
   10   continue
      end if
c
      if (etasw .ge. 2 .and. eta.ne. 0.0) then
        jcrit1=jcrit0*delbz
        jcrit2=jcrit1**2
        eta0 = eta*(exp(-time/20.0)-exp(-time/2.0))
c        eta0 = 0.0
        eta1 = 0.000
        do 50 iz = 1,nzz
	do 50 iy = 1,ny
	do 50 ix = 1,nx
	  hilfy(ix,iy,iz) = res(ix,iy,iz)
	  res(ix,iy,iz) = 0.
	  hilfx(ix,iy,iz) = 0.
   50   continue
c hilf=j**2 calculated in lax and leap

c    model 1
        if (etasw .eq. 2) then
          do 60 iz = 2,nzz-1
          do 60 iy = 2,ny1
	  ixanf = 2 + mod(iz+iy+igrid,2)
	  do 60 ix = ixanf,nx1,2
            if (hilf(ix,iy,iz) .gt. jcrit1)  
     +          res(ix,iy,iz) = eta*sqrt(hilf(ix,iy,iz)-jcrit1)
   60     continue
        endif
        if (etasw .ge. 3) then
          do 70 iz = 2,nzz-1
          do 70 iy = 2,ny1
	  ixanf = 2 + mod(iz+iy+igrid,2)
	  do 70 ix = ixanf,nx1,2
            if (hilf(ix,iy,iz) .gt. jcrit2)  
     +          res(ix,iy,iz) = eta*(hilf(ix,iy,iz)-jcrit2)
   70     continue
        endif

        call resbd

        do 80 iz = 2,nzz-1
        do 80 iy = 2,ny1
	ixanf = 3 - mod(iz+iy+igrid,2)
	do 80 ix = ixanf,nx1,2
	  res(ix,iy,iz)  = 0.33333*( (meanpx(ix)*res(ix+1,iy,iz)
     +                        + meanmx(ix)*res(ix-1,iy,iz) )
     +                   + (meanpy(iy)*res(ix,iy+1,iz)
     +                        + meanmy(iy)*res(ix,iy-1,iz) )
     +                   + (meanpz(iz)*res(ix,iy,iz+1)
     +                        + meanmz(iz)*res(ix,iy,iz-1) ) )
   80   continue
        do 85 iz = 1,nzz
	do 85 iy = 1,ny
	do 85 ix = 1,nx
           res(ix,iy,iz) = res(ix,iy,iz) + eta0*prof(ix,iy,iz) + eta1
   85   continue

        call resbd

	do 89 ix = 1,nx
	  ppp1(ix) =  (1.-0.75/cosh(x(ix)/2.))
          ppp2(ix) =  0.25/cosh(x(ix)/2.)
   89   continue
        do 90 iz = 2,nzz-1
	do 90 iy = 2,ny1
	do 90 ix = 2,nx1
	  hilfx(ix,iy,iz) =  ppp1(ix)*res(ix,iy,iz)
     +              + ppp2(ix)*(meanpx(ix)*res(ix+1,iy,iz)
     +                        + meanmx(ix)*res(ix-1,iy,iz)
     +                        + meanpy(iy)*res(ix,iy+1,iz)
     +                        + meanmy(iy)*res(ix,iy-1,iz)
     +                        + meanpz(iz)*res(ix,iy,iz+1)
     +                        + meanmz(iz)*res(ix,iy,iz-1) )
   90   continue
        epar1=1.0-dt
        do 100 iz = 2,nzz-1
	do 100 iy = 2,ny1
	do 100 ix = 2,nx1
	    res(ix,iy,iz) = epar1*hilfy(ix,iy,iz)+dt*hilfx(ix,iy,iz)
  100   continue
c
        call resbd
      end if
c
      return
      end
	subroutine binin
c***********************************************************************
	include 'misflin'
c
	integer           ier,anzx,anzy,anzz,anzpos,anzrow
        character*26     comagt
c................................................................
         comagt='-a magtap   -F f77 -N ieee'
c  uncomment 'call asnunit' for cray	 
         if (ieein) then 
           call asnunit(14,comagt,ier)
         else
	   open (14,file='magtap',form='unformatted')
	 end if
	 read(14) anzx,anzy,anzz,time,anzpos,anzrow
	 if ( anzx.ne.nx .or. anzy.ne.ny .or. anzz.ne.nz
     +        .or. anzpos.ne.npos .or. anzrow.ne.nrow )  then
	     write(26,11)  anzx,anzy,anzz,anzpos,anzrow
	     fini = .true.
	     return
         end if
         read(14)   x,difx,ddifx,ddifpx,ddifmx,meanpx,meanmx,
     +              y,dify,ddify,ddifpy,ddifmy,meanpy,meanmy,
     +              z,difz,ddifz,ddifpz,ddifmz,meanpz,meanmz
         read(14)   bx,by,bz
         read(14)   sx,sy,sz
         read(14)   rho,u,res,prof
         read(14)   istep,isafe,iend,intart,ivisrho,nsmooth,
     +              iferror,igrid,idiag,nzz,cval,
     +              ferror,dt,gamma,eta,visx,visy,visz,
     +              nsim,bsim,vsim,psim,lsim,tsim,
     +              psi,phi,bmsh,pmsp,kappa,delrho,xrho,dxrho,
     +              delbz,bzmsp,bzmsh,pmsh,betasp,betash,
     +              vamsp,vamsh,csmsp,csmsh,
     +              perio,lsym,kbx,kby,kbz,ksx,ksy,ksz,
     +              arho,au,abx,aby,abz,asx,asy,asz,
     +              xmin,xmax,ymin,ymax,zmin,zmax,aequi,zentr,eps,
     +              rho0,p0,u0,vx0,vy0,vz0,bx0,by0,bz0,
     +              rho1,p1,u1,vx1,vy1,vz1,bx1,by1,bz1,
     +              rhoprof,uprof,
     +              bxprof,byprof,bzprof,sxprof,syprof,szprof,
     +              rhobd,pbd,ubd,
     +              vxbd,vybd,vzbd,bxbd,bybd,bzbd,
     +              tsatout,isat,timesat,
     +              xsat,ysat,zsat,vxsat,vysat,vzsat
c  uncomment 'call assig'  for cray	 
	  if (ieein)  call assign('assign -R')
	 close (14)
c
   11 format(1x,' !!warning!!!',
     +       /1x,'Grid parameters: nx=',i4,'  ny=',i4,
     +           '  ny=',i4,
     +       /1x,'or satparameters: npos=',i4,'  nrow=',i4,
     +       /1x,'not consistent with program parameters',
     +       /1x,'!!!warning!!!')

      return
      end
	subroutine binout(ind)
c***********************************************************************
	include 'misflin'
c
        integer          ind,ier
	character*8      magtapn
        character*26     comagt
c................................................................
         write(magtapn,1969) ind
         write(comagt,1970) ind
 1969    format('magtap',i2.2)
 1970    format('-a magtap',i2.2,' -F f77 -N ieee')
c  uncomment 'call asnunit' for cray	 
         if (ieeout) then 
           call asnunit(14,comagt,ier)
         else
	   open (14,file=magtapn,form='unformatted')
	 end if
         write(14)   nx,ny,nz,time,npos,nrow
         write(14)   x,difx,ddifx,ddifpx,ddifmx,meanpx,meanmx,
     +               y,dify,ddify,ddifpy,ddifmy,meanpy,meanmy,
     +               z,difz,ddifz,ddifpz,ddifmz,meanpz,meanmz
         write(14)   bx,by,bz
         write(14)   sx,sy,sz
         write(14)   rho,u,res,prof
         write(14)   istep,isafe,iend,intart,ivisrho,nsmooth,
     +              iferror,igrid,idiag,nzz,cval,
     +              ferror,dt,gamma,eta,visx,visy,visz,
     +              nsim,bsim,vsim,psim,lsim,tsim,
     +              psi,phi,bmsh,pmsp,kappa,delrho,xrho,dxrho,
     +              delbz,bzmsp,bzmsh,pmsh,betasp,betash,
     +              vamsp,vamsh,csmsp,csmsh,
     +              perio,lsym,kbx,kby,kbz,ksx,ksy,ksz,
     +              arho,au,abx,aby,abz,asx,asy,asz,
     +              xmin,xmax,ymin,ymax,zmin,zmax,aequi,zentr,eps,
     +              rho0,p0,u0,vx0,vy0,vz0,bx0,by0,bz0,
     +              rho1,p1,u1,vx1,vy1,vz1,bx1,by1,bz1,
     +              rhoprof,uprof,
     +              bxprof,byprof,bzprof,sxprof,syprof,szprof,
     +              rhobd,pbd,ubd,
     +              vxbd,vybd,vzbd,bxbd,bybd,bzbd,
     +              tsatout,isat,timesat,
     +              xsat,ysat,zsat,vxsat,vysat,vzsat
c  uncomment following line for cray	 
	 if (ieeout)  call assign('assign -R')
	 close (14)
	 
      return
      end
      subroutine initflux
c***********************************************************************
	include 'misflin'
c.................................................................
        open (16,file='flpar')
        read (16,1) wsmooth,mhdaus,baus,jaus,binaus,mirror
        read (16,3) xmineu,xmaxneu,ymineu,ymaxneu,zmineu,zmaxneu
        read (16,1) unif(1),unif(2),unif(3)
        read (16,1) zenew(1),zenew(2),zenew(3)
        read (16,3) epsnew(1),epsnew(2),epsnew(3)
        write(*,3) epsnew(1),epsnew(2),epsnew(3)
        read (16,1) mpause
        write(*,1) mpause
        read (16,3) deltfl,xfmin,xfmax,yfmin,yfmax,zfmin,zfmax
        write(*,3) xfmin,xfmax,yfmin,yfmax,zfmin,zfmax
        read (16,2) nplane,anzint
        write(*,2) nplane,anzint
        read (16,3) dels
        write(*,3) dels
        if (nplane.gt.5) nplane=1
        do 5 it = 1, nplane
          read (16,3) start00(1,it),start00(2,it),start00(3,it)
          read (16,3) start01(1,it),start01(2,it),start01(3,it)
          read (16,3) start10(1,it),start10(2,it),start10(3,it)
    5   continue
	call gridh
        write(6,*) 'in init vor parfl'
	call parfl
	call parflux
	write(6,*) 'in init hinter parfl'
      close (16)
    1 format(1x,l10)
    2 format(1x,i9)
    3 format(1x,f14.5)
c
      return
      end
      subroutine orgflux
c***********************************************************************
	include 'misflin'
c
       character*9      mflcray
       character*27     mfluxn
c.................................................................
       do 100 iff = 1,nplane
          if (ieeout .eq. 1) then
           write(mfluxn,1) iff,nflux
           call asnunit(22,mfluxn,ier)
          else
           write(mflcray,2) iff,nflux
	   open (22,file=mflcray,form='unformatted')
	  endif
	 
	  write(6,*) 'in org vor fluxfun'
	  call fluxfun(iff)
	  if (ieeout .eq. 1)   call assign('assign -R')
          close (22)
  100  continue
       nflux=nflux+1
        
    1 format('-a mflux',i1.1,i3.3,' -F f77 -N ieee')
    2 format('mflux',i1.1,i3.3)
c
      return
      end
      subroutine parflux
c***********************************************************************
	include 'misflin'
c
      include 'misflin'
c   ------------------------------------
c    ausgabenkram
c   ------------------------------------
      open (21,file='magparflux')
      write(21,1000)  istep,time
      write(21,1001)  nx,ny,nz
c
      write(21,1011)
      do 101 ix = 1,nx
	write(21,1021) x(ix),difx(ix),ddifx(ix),ddifpx(ix),ddifmx(ix),
     +                 meanpx(ix),meanmx(ix)
  101 continue
      write(21,1012)
      do 102 iy = 1,ny
	write(21,1021) y(iy),dify(iy),ddify(iy),ddifpy(iy),ddifmy(iy),
     +                 meanpy(iy),meanmy(iy)
  102 continue
      write(21,1013)
      do 103 iz = 1,nz
	write(21,1021) z(iz),difz(iz),ddifz(iz),ddifpz(iz),ddifmz(iz),
     +                 meanpz(iz),meanmz(iz)
  103 continue
c
      write(21,1081) 
      write(21,1082) wsmooth, fluxcal
      write(21,1083) mhdaus,baus,jaus,mirror
      write(21,1084) xmineu,ymineu,zmineu,xmaxneu,ymaxneu,zmaxneu
      write(21,1085)  unif(1),unif(2),unif(3),
     +              zenew(1),zenew(2),zenew(3),
     +              epsnew(1),epsnew(2),epsnew(3)
c
c
      if (fluxcal) then
        write(21,2099)
        write (21,2101) xfmin,xfmax
        write (21,2102) yfmin,yfmax
        write (21,2103) zfmin,zfmax
        write (21,2104) nxf,nyf,nzf
        write (21,2105) anzint,dels
        write (21,2106) nplane
        do 501 ip = 1,nplane
          write (21,2107) ip,start00(1,ip),start00(2,ip),start00(3,ip)
          write (21,2108) start01(1,ip),start01(2,ip),start01(3,ip)
          write (21,2108) start10(1,ip),start10(2,ip),start10(3,ip)
  501   continue
        write (21,2109) 
        do 550 ix = 1,nxf
	  write(21,2201) ix,xf(ix)
  550   continue
        do 560 iy = 1,nyf
	  write(21,2202) iy,yf(iy)
  560   continue
        do 570 iz = 1,nzf
	  write(21,2203) iz,zf(iz)
  570   continue
      end if
      
      close (21)
c
 1000 format(/1x,'program parameter(intstep=',i5,',  time=', 
     +            f7.2,'): ',
     +       /1x,'-------------------------------------------------')
c
 1001 format(/1x,'nx     =',i5,3x,'ny      =',i5,3x,'nz     =',i5,3x)
c
 1011 format(/4x,'  x         dif       ddif      ddifp     ddifm',
     +           '     meanp     meanm ')
 1012 format(/4x,'  y         dif       ddif      ddifp     ddifm',
     +           '     meanp     meanm ')
 1013 format(/4x,'  z         dif       ddif      ddifp     ddifm',
     +           '     meanp     meanm ')
 1021 format(1x,7f10.3)
c
 1081 format (/1x, '************************',
     +             'Parameter of New Field',
     +              '***********************')
 1082 format (1x, 'With smoothing:  ',l3,
     +        1x, '  Calculate flux surface along field lines: ',l3)
 1083 format (1x, 'Output MHD: ',l3, '  Output B field: ',l3,
     +             '  Output J field: ',l3,'  Mirror at z=0: ',l3)
 1084 format(1x,'xminew =',f9.3,3x,'yminew =',f9.3,3x,'zminew =',f9.3,
     +      /1x,'xmaxew =',f9.3,3x,'ymaxew =',f9.3,3x,'zmaxew =',f9.3)
 1085 format(/1x,'uniform grid in       x:',l3,8x,'y:',l3,
     +        8x,'z:',l3,
     +       /1x,'max resolution in   x:',l3,
     +        8x,'y:',l3,8x,'z:',l3,
     +       /1x,'with grid spacing:         x:',f7.3,
     +        4x,'y:',f7.3,4x,'z:',f7.3)
 1101 format('  xneu(',i3,') = ',f8.4,'  dxneu(',i3,') = ',f8.4)
 1102 format('  yneu(',i3,') = ',f8.4,'  dyneu(',i3,') = ',f8.4)
 1103 format('  zneu(',i3,') = ',f8.4,'  dzneu(',i3,') = ',f8.4)
c
 2099 format (/1x, '*******************',
     +             'Flux Surface Computation',
     +             '*******************')
 2101 format (1x,'boundaries in x: xmin =',f10.6,' xmax =',f10.6)
 2102 format (1x,'boundaries in y: ymin =',f10.6,' ymax =',f10.6)
 2103 format (1x,'boundaries in z: zmin =',f10.6,' zmax =',f10.6)
 2104 format (1x,'no. of gridpoints in x =',i4,
     +           '  in y:',i4,'  in z:',i4)
 2105 format (1x,'no. of field integrations:',i5,'  int step:',f8.4)
 2106 format (1x,'no. datasets',i3)
 2107 format (1x,'corners for startlocation (',i2,'): x=',f9.2,
     +           '  y=',f9.2,'  y=',f9.2)
 2108 format (1x,'                                x=',f9.2,
     +           '  y=',f9.2,'  y=',f9.2)
 2109 format (1x,'Grid Coordinates:')
 2201 format('  xf(',i3,') = ',f8.4)
 2202 format('  yf(',i3,') = ',f8.4)
 2203 format('  zf(',i3,') = ',f8.4)
c
      return
      end
      subroutine gridh
c***********************************************************************
	include 'misflin'
c
      include 'misflin'
c.................................................................
      do 100 ix = 1, nx-1
	x3d(ix)  = x(ix)
	dx3d(ix) = x(ix+1)-x(ix)
  100 continue
      x3d(nx)  = x(nx)
      do 200 iy = 1, ny-1
        y3d(iy) = y(iy)
	dy3d(iy) = y(iy+1)-y(iy)
  200 continue
      y3d(ny)  = y(ny)
      do 300 iz = 1, nz-1
	z3d(iz) = z(iz)
	dz3d(iz) = z(iz+1)-z(iz)
  300 continue
      z3d(nz)  = z(nz)
c
      return
      end
      subroutine parfl
c***********************************************************************
	include 'misflin'
c
      include 'misflin'
c.................................................................
      xfsize = xfmax-xfmin
      yfsize = yfmax-yfmin
      zfsize = zfmax-zfmin
      delxf = xfsize/(nxf-1)
      delyf = yfsize/(nyf-1)
      delzf = zfsize/(nzf-1)
      if (anzint*dels .lt. 3.*zmax) write(*,11) dels, anzint, zmax
      if (anzint*dels .lt. 2.*ymax) write(*,12) dels, anzint, ymax
      do 100 ix = 2, nxf-1
	xf(ix) = xfmin+delxf*(ix-1)
  100 continue
      xf(1)   = xfmin
      xf(nxf) = xfmax
      do 200 iy = 2, nyf-1
	yf(iy) = yfmin+delyf*(iy-1)
  200 continue
      yf(1)   = yfmin
      yf(nyf) = yfmax
      do 300 iz = 2, nzf-1
	zf(iz) = zfmin+delzf*(iz-1)
  300 continue
      zf(1)   = zfmin
      zf(nzf) = zfmax
c
   11 format(/1x,'!!!warning!!!',
     +       /1x,'integration parameter: dels (',f10.4,
     +           ')  times anz intsteps ( ',i4,')  ',
     +       /1x,'is less than 3 times zmax (',f10.4,')')
   12 format(/1x,'!!!warning!!!',
     +       /1x,'integration parameter: dels (',f10.4,
     +           ')  times anz intsteps ( ',i4,')  ',
     +       /1x,'is less than twice ymax (',f10.4,')')
      return
      end
      subroutine fluxfun(iff)
c***********************************************************************
	include 'misflin'
c
      integer  iff,it,ic,ix,iy,iz,ifl1,ifl2,isumtot
      real     dst1(3),dst2(3),strt(3)
c...............................................
      isumtot=0
      do 2 ix = 1, nxf
       do 2 iy = 1, nyf
        do 2 iz = 1, nzf
    2 ff(ix,iy,iz)=0.0
      do 4 indf = 1, nindf
        iw(indf)=0
        fw(indf)=0.0
    4 continue
c
c        write(*,*) 'in fluxfun'
      do 10 ic = 1, 3
         dst1(ic) = start01(ic,iff)/(nfl1-1)
         dst2(ic) = start10(ic,iff)/(nfl2-1)
   10 continue
      do 100 ifl2 = 1, nfl2
	  strt(1) = start00(1,iff)+dst2(1)*(ifl2-1)
	  strt(2) = start00(2,iff)+dst2(2)*(ifl2-1)
	  strt(3) = start00(3,iff)+dst2(3)*(ifl2-1)
        do 20 ifl1 = 1, nfl1
	  rx(ifl1) = strt(1)+dst1(1)*(ifl1-1)
	  ry(ifl1) = strt(2)+dst1(2)*(ifl1-1)
	  rz(ifl1) = strt(3)+dst1(3)*(ifl1-1)
   20   continue
        inatot=.true.
        do 30 ii = 1, nfl1
   30     inarea(ii) = 1
        delint=dels
        lastst=1
	do 60 it = 1, anzint
          if (mpause.and.inatot) call bintmp
          if (.not.mpause.and.inatot) call bintarc
          if (inatot) lastst=it
          isumtot=isumtot+isum
   60   continue
c        if (inatot) then
c          write(*,11) iff, ifl2
c          do 70 ii = 1, nfl1
c   70     if (inarea(ii) .eq. 1) write(*,12) ii,rx(ii),ry(ii),rz(ii)
c        endif
c        write(*,13) ifl2, lastst, isumtot

        do 120 ifl1 = 1, nfl1
	  rx(ifl1) = strt(1)+dst1(1)*(ifl1-1)
	  ry(ifl1) = strt(2)+dst1(2)*(ifl1-1)
	  rz(ifl1) = strt(3)+dst1(3)*(ifl1-1)
  120   continue
        inatot=.true.
        do 130 ii = 1, nfl1
  130     inarea(ii) = 1
        delint=-dels
        lastst=1
	do 160 it = 1, anzint
          if (mpause.and.inatot) call bintmp
          if (.not.mpause.and.inatot) call bintarc
          if (inatot) lastst=it
          isumtot=isumtot+isum
  160   continue
c        if (inatot) then
c          write(*,11) iff, ifl2
c          do 170 ii = 1, nfl1
c  170     if (inarea(ii) .eq. 1) write(*,12) ii,rx(ii),ry(ii),rz(ii)
c        endif
c        write(*,13) ifl2, lastst, isumtot

  100 continue
c      write(*,13) ifl2, lastst, isumtot
c            
      call wflux(iaa,faa)
c
   11 format(/1x,'!!!warning!!!',
     +           ' Integration terminated in area with',
     +           ' following parameters:',
     +       /1x,'   Plane (',i3,')  ifl2 =',i4)
   12 format(1x,'     ifl1 =',i4,'  x=',f8.2,'  y=',f8.2,'  z=',f8.2)
   13 format(1x,'SET ifl2 =',i4,' last step:',i4,' total isteps',i8)
   14 format(1x,'NUMBER OF NONZERO ELEMENTS IN FF:' i8)
c
      return
      end
      subroutine wflux
c***********************************************************************
	include 'misflin'
c
      integer   in,ix,iy,iz,indf
c.................................................................
      write(22) nxf,nyf,nzf,time,nfl1,nfl2
      write(22) xf,yf,zf

      do 100 in = 1, nindf
        iw(in)=0
        fw(in)=0.0
  100 continue
      
      indf=0
      do 200 ix = 1, nxf
       do 200 iy = 1, nyf
        do 200 iz = 1, nzf
         if (ff(ix,iy,iz).ge.0.25)   indf=indf+1
  200 continue
      nonon0=indf
c      write(*,*) 'No. non zero>0.25 nindf', nonon0
      indf=0
      do 210 ix = 1, nxf
       do 210 iy = 1, nyf
        do 210 iz = 1, nzf
         if (ff(ix,iy,iz).ge.0.25 .and. indf.le.nindf-1) then
           indf=indf+1
           iw(indf)=ix+(iy-1)*nxf+(iz-1)*nxf*nyf
           fw(indf)=ff(ix,iy,iz)
c           if (indf.lt.11 .or. (nonon0-indf).lt.11) 
c     +         write(*,*) indf,ix,iy,iz,iw(indf),fw(indf)
         endif      
  210 continue
      if (nonon0.gt.nindf) 
     +    write(*,*) 'No. non zero > nindf', nonon0,nindf
      if (nonon0.gt.nindf) nonon0=nindf
      write(22) nonon0,nindf
      write(22) iw,fw
c
      return
      end
      subroutine bintmp 
c***********************************************************************
	include 'misflin'
c
      integer   ii, indfast, zpmnull(nfl1), inarf(nfl1)
      integer   iz0,iz1,iv00,iv01,iv10,iv11
      integer   ix(nfl1), iy(nfl1), iz(nfl1), iv(nfl1),
     +          ixf(nfl1), iyf(nfl1), izf(nfl1)
      real      hxu(nfl1), hyu(nfl1), hzu(nfl1),
     +          hxo(nfl1), hyo(nfl1), hzo(nfl1)
      real      hxf(nfl1), hyf(nfl1), hzf(nfl1),
     +          hxf1(nfl1), hyf1(nfl1), hzf1(nfl1)
      real      ttx(nfl1), tty(nfl1), ttz(nfl1)
      real      ttt(nfl1)
      real      rxa(nfl1), rya(nfl1), rza(nfl1)
c
c  Integration of field lines using a 2 step 
c           lax wendroff scheme.
c   >indfast< determines grid indices for field interpolation.
c   Care is required for special symmetries!
c   Symmetries considered are 
c   (a) linesymmetry along x axis for mpause (mpause=true)
c   (b) linesymmetry along y axis for arc (mpause=false)
c   Use of symmetry requires mirror=true
c.........................................
c I. Step
c
      do 20 ii = 1, nfl1
	if (rx(ii).ge.x3d(nx-1))   inarea(ii) = 0
	if (rx(ii).le.x3d(   2))   inarea(ii) = 0
	if (ry(ii).ge.y3d(ny-1))   inarea(ii) = 0
	if (ry(ii).le.y3d(   2))   inarea(ii) = 0
	if (rz(ii).ge.z3d(nz-1))   inarea(ii) = 0
	if (mirror .and. rz(ii).le.-z3d(nz-1))  inarea(ii) = 0
	if (.not.mirror .and. rz(ii).le.z3d(2))  inarea(ii) = 0
   20 continue
      isum=0
      do 22 ii = 1, nfl1
        isum=isum+inarea(ii)
   22 continue
      if (isum.eq.0) then
        inatot=.false.
        return
      endif

      do 30 ii = 1, nfl1
	if (rz(ii) .lt. 0) then
	  zpmnull(ii) = -1
	else
	  zpmnull(ii) = 1
	end if
   30 continue

      do 50 ii = 1, nfl1
       if (inarea(ii).eq.1) then
	ix(ii) = indfast (rx(ii),x3d,nx)
	iy(ii) = indfast (zpmnull(ii)*ry(ii),y3d,ny)
        iz(ii) = indfast (zpmnull(ii)*rz(ii),z3d,nz)
        iv(ii) = ix(ii) + (iy(ii)-1)*nx
	hxo(ii) = (rx(ii)-x3d(ix(ii))) / dx3d(ix(ii))
        hyo(ii) = (zpmnull(ii)*ry(ii)-y3d(iy(ii))) / dy3d(iy(ii))
        hzo(ii) = (zpmnull(ii)*rz(ii)-z3d(iz(ii))) / dz3d(iz(ii))
c        if (iy(ii).eq.ny-1) write(*,*) 
c     +                'iy=ny-1',ii,ry(ii),hyo(ii),y3d(iy(ii))
       endif
   50 continue
      do 60 ii = 1, nfl1
       if (inarea(ii).eq.1) then
	 hxu(ii) = 1.0 - hxo(ii)
	 hyu(ii) = 1.0 - hyo(ii)
	 hzu(ii) = 1.0 - hzo(ii)
       endif
   60 continue
c
      do 90 ii = 1, nfl1
       if (inarea(ii).eq.1) then
	 ix0 = ix(ii)
	 ix1 = ix(ii) + 1
	 iy0 = iy(ii)
	 iy1 = iy(ii) + 1
	 iz0 = iz(ii)
	 iz1 = iz(ii) + 1
	 ttx(ii) = zpmnull(ii)*(
     +               hxu(ii) * hyu(ii) * hzu(ii) * bx(ix0,iy0,iz0)
     +             + hxu(ii) * hyu(ii) * hzo(ii) * bx(ix0,iy0,iz1)
     +             + hxu(ii) * hyo(ii) * hzu(ii) * bx(ix0,iy1,iz0)
     +             + hxu(ii) * hyo(ii) * hzo(ii) * bx(ix0,iy1,iz1)
     +             + hxo(ii) * hyu(ii) * hzu(ii) * bx(ix1,iy0,iz0)
     +             + hxo(ii) * hyu(ii) * hzo(ii) * bx(ix1,iy0,iz1)
     +             + hxo(ii) * hyo(ii) * hzu(ii) * bx(ix1,iy1,iz0)
     +             + hxo(ii) * hyo(ii) * hzo(ii) * bx(ix1,iy1,iz1) )
	 tty(ii) = ( hxu(ii) * hyu(ii) * hzu(ii) * by(ix0,iy0,iz0)
     +             + hxu(ii) * hyu(ii) * hzo(ii) * by(ix0,iy0,iz1)
     +             + hxu(ii) * hyo(ii) * hzu(ii) * by(ix0,iy1,iz0)
     +             + hxu(ii) * hyo(ii) * hzo(ii) * by(ix0,iy1,iz1)
     +             + hxo(ii) * hyu(ii) * hzu(ii) * by(ix1,iy0,iz0)
     +             + hxo(ii) * hyu(ii) * hzo(ii) * by(ix1,iy0,iz1)
     +             + hxo(ii) * hyo(ii) * hzu(ii) * by(ix1,iy1,iz0)
     +             + hxo(ii) * hyo(ii) * hzo(ii) * by(ix1,iy1,iz1) )
	 ttz(ii) = ( hxu(ii) * hyu(ii) * hzu(ii) * bz(ix0,iy0,iz0)
     +             + hxu(ii) * hyu(ii) * hzo(ii) * bz(ix0,iy0,iz1)
     +             + hxu(ii) * hyo(ii) * hzu(ii) * bz(ix0,iy1,iz0)
     +             + hxu(ii) * hyo(ii) * hzo(ii) * bz(ix0,iy1,iz1)
     +             + hxo(ii) * hyu(ii) * hzu(ii) * bz(ix1,iy0,iz0)
     +             + hxo(ii) * hyu(ii) * hzo(ii) * bz(ix1,iy0,iz1)
     +             + hxo(ii) * hyo(ii) * hzu(ii) * bz(ix1,iy1,iz0)
     +             + hxo(ii) * hyo(ii) * hzo(ii) * bz(ix1,iy1,iz1) )
       endif
  90  continue
c
      do 95 ii = 1, nfl1
       if (inarea(ii).eq.1) then
         ttt(ii)=sqrt(ttx(ii)*ttx(ii)+tty(ii)*tty(ii)+ttz(ii)*ttz(ii))
       endif
  95  continue
      do 100 ii = 1, nfl1
       if (inarea(ii).eq.1) then
	rxa(ii) = rx(ii) + 0.5*delint*ttx(ii)/ttt(ii)
	rya(ii) = ry(ii) + 0.5*delint*tty(ii)/ttt(ii)
	rza(ii) = rz(ii) + 0.5*delint*ttz(ii)/ttt(ii)
c	 if (tty(ii).le.0.1) write(*,*) 'tty<0.1', ii,tty(ii)
       endif
  100 continue
c
c II. Step
c
      do 120 ii = 1, nfl1
       if (inarea(ii).eq.1) then
	if (rxa(ii).ge.x3d(nx-1))   inarea(ii) = 0
	if (rxa(ii).le.x3d(   2))   inarea(ii) = 0
	if (rya(ii).ge.y3d(ny-1))   inarea(ii) = 0
	if (rya(ii).le.y3d(   2))   inarea(ii) = 0
	if (rza(ii).ge.z3d(nz-1))   inarea(ii) = 0
	if (mirror .and. rza(ii).le.-z3d(nz-1))  inarea(ii) = 0
	if (.not.mirror .and. rza(ii).le.z3d(2))  inarea(ii) = 0
       endif
  120 continue
      isum=0
      do 122 ii = 1, nfl1
        isum=isum+inarea(ii)
  122 continue
      if (isum.eq.0) then
        inatot=.false.
        return
      endif
      do 130 ii = 1, nfl1
	 if (rza(ii) .lt. 0) then
	    zpmnull(ii) = -1
	 else
	    zpmnull(ii) = 1
	 end if
  130 continue
c
c   Scalar flux field:
c
      do 103 ii = 1, nfl1
  103   inarf(ii) = 1
      do 105 ii = 1, nfl1
	if (rxa(ii).ge.xfmax .or. rxa(ii).le.xfmin) inarf(ii) = 0
	if (rya(ii).ge.yfmax .or. rya(ii).le.yfmin) inarf(ii) = 0
	if (rza(ii).ge.zfmax .or. rza(ii).le.zfmin) inarf(ii) = 0
  105 continue
      do 110 ii = 1, nfl1
       if (inarf(ii).eq.1) then
         ixf(ii)  = int((rxa(ii)-xfmin)/xfsize*(nxf-1))+1
         iyf(ii)  = int((rya(ii)-yfmin)/yfsize*(nyf-1))+1
         izf(ii)  = int((rza(ii)-zfmin)/zfsize*(nzf-1))+1
         hxf(ii)  = (rxa(ii)-xf(ixf(ii)))/delxf
         hyf(ii)  = (rya(ii)-yf(iyf(ii)))/delyf
         hzf(ii)  = (rza(ii)-zf(izf(ii)))/delzf
         if (hxf(ii).lt.0) hxf(ii)=0.0
         if (hyf(ii).lt.0) hyf(ii)=0.0
         if (hzf(ii).lt.0) hzf(ii)=0.0
         hxf1(ii) = 1.-hxf(ii)
         hyf1(ii) = 1.-hyf(ii)
         hzf1(ii) = 1.-hzf(ii)
       endif
  110 continue
      do 115 ii = 1, nfl1
       if (inarf(ii).eq.1) then
         ff(ixf(ii),iyf(ii),izf(ii))=
     +     ff(ixf(ii),iyf(ii),izf(ii))+hxf1(ii)*hyf1(ii)*hzf1(ii)
         ff(ixf(ii)+1,iyf(ii),izf(ii))=
     +     ff(ixf(ii)+1,iyf(ii),izf(ii))+hxf(ii)*hyf1(ii)*hzf1(ii)
         ff(ixf(ii),iyf(ii)+1,izf(ii))=
     +     ff(ixf(ii),iyf(ii)+1,izf(ii))+hxf1(ii)*hyf(ii)*hzf1(ii)
         ff(ixf(ii),iyf(ii),izf(ii)+1)=
     +     ff(ixf(ii),iyf(ii),izf(ii)+1)+hxf1(ii)*hyf1(ii)*hzf(ii)
         ff(ixf(ii)+1,iyf(ii)+1,izf(ii))=
     +     ff(ixf(ii)+1,iyf(ii)+1,izf(ii))+hxf(ii)*hyf(ii)*hzf1(ii)
         ff(ixf(ii),iyf(ii)+1,izf(ii)+1)=
     +     ff(ixf(ii),iyf(ii)+1,izf(ii)+1)+hxf1(ii)*hyf(ii)*hzf(ii)
         ff(ixf(ii)+1,iyf(ii),izf(ii)+1)=
     +     ff(ixf(ii)+1,iyf(ii),izf(ii)+1)+hxf(ii)*hyf1(ii)*hzf(ii)
         ff(ixf(ii)+1,iyf(ii)+1,izf(ii)+1)=
     +     ff(ixf(ii)+1,iyf(ii)+1,izf(ii)+1)+hxf(ii)*hyf(ii)*hzf(ii)
       endif
  115 continue
c
      do 150 ii = 1, nfl1
       if (inarea(ii).eq.1) then
	 ix(ii) = indfast (rxa(ii),x3d,nx)
	 iy(ii) = indfast (zpmnull(ii)*rya(ii),y3d,ny)
	 iz(ii) = indfast (zpmnull(ii)*rza(ii),z3d,nz)
	 iv(ii) = ix(ii) + (iy(ii)-1)*nx
	 hxo(ii) = (rxa(ii)-x3d(ix(ii))) / dx3d(ix(ii))
	 hyo(ii) = (zpmnull(ii)*rya(ii)-y3d(iy(ii))) / dy3d(iy(ii))
	 hzo(ii) = (zpmnull(ii)*rza(ii)-z3d(iz(ii))) / dz3d(iz(ii))
c        if (iy(ii).eq.ny-1) write(*,*) 
c     +                'iy=ny-1',ii,ry(ii),hyo(ii),y3d(iy(ii))
       endif
 150  continue
      do 160 ii = 1, nfl1
       if (inarea(ii).eq.1) then
	 hxu(ii) = 1.0 - hxo(ii)
	 hyu(ii) = 1.0 - hyo(ii)
	 hzu(ii) = 1.0 - hzo(ii)
       endif
 160  continue
c
      do 190 ii = 1, nfl1
       if (inarea(ii).eq.1) then
	 ix0 = ix(ii)
	 ix1 = ix(ii) + 1
	 iy0 = iy(ii)
	 iy1 = iy(ii) + 1
	 iz0 = iz(ii)
	 iz1 = iz(ii) + 1
	 ttx(ii) = zpmnull(ii)*(
     +               hxu(ii) * hyu(ii) * hzu(ii) * bx(ix0,iy0,iz0)
     +             + hxu(ii) * hyu(ii) * hzo(ii) * bx(ix0,iy0,iz1)
     +             + hxu(ii) * hyo(ii) * hzu(ii) * bx(ix0,iy1,iz0)
     +             + hxu(ii) * hyo(ii) * hzo(ii) * bx(ix0,iy1,iz1)
     +             + hxo(ii) * hyu(ii) * hzu(ii) * bx(ix1,iy0,iz0)
     +             + hxo(ii) * hyu(ii) * hzo(ii) * bx(ix1,iy0,iz1)
     +             + hxo(ii) * hyo(ii) * hzu(ii) * bx(ix1,iy1,iz0)
     +             + hxo(ii) * hyo(ii) * hzo(ii) * bx(ix1,iy1,iz1) )
	 tty(ii) = ( hxu(ii) * hyu(ii) * hzu(ii) * by(ix0,iy0,iz0)
     +             + hxu(ii) * hyu(ii) * hzo(ii) * by(ix0,iy0,iz1)
     +             + hxu(ii) * hyo(ii) * hzu(ii) * by(ix0,iy1,iz0)
     +             + hxu(ii) * hyo(ii) * hzo(ii) * by(ix0,iy1,iz1)
     +             + hxo(ii) * hyu(ii) * hzu(ii) * by(ix1,iy0,iz0)
     +             + hxo(ii) * hyu(ii) * hzo(ii) * by(ix1,iy0,iz1)
     +             + hxo(ii) * hyo(ii) * hzu(ii) * by(ix1,iy1,iz0)
     +             + hxo(ii) * hyo(ii) * hzo(ii) * by(ix1,iy1,iz1) )
	 ttz(ii) = ( hxu(ii) * hyu(ii) * hzu(ii) * bz(ix0,iy0,iz0)
     +             + hxu(ii) * hyu(ii) * hzo(ii) * bz(ix0,iy0,iz1)
     +             + hxu(ii) * hyo(ii) * hzu(ii) * bz(ix0,iy1,iz0)
     +             + hxu(ii) * hyo(ii) * hzo(ii) * bz(ix0,iy1,iz1)
     +             + hxo(ii) * hyu(ii) * hzu(ii) * bz(ix1,iy0,iz0)
     +             + hxo(ii) * hyu(ii) * hzo(ii) * bz(ix1,iy0,iz1)
     +             + hxo(ii) * hyo(ii) * hzu(ii) * bz(ix1,iy1,iz0)
     +             + hxo(ii) * hyo(ii) * hzo(ii) * bz(ix1,iy1,iz1) )
       endif
 190  continue
c
      do 195 ii = 1, nfl1
       if (inarea(ii).eq.1) then
         ttt(ii)=sqrt(ttx(ii)*ttx(ii)+tty(ii)*tty(ii)+ttz(ii)*ttz(ii))
       endif
 195  continue
      do 200 ii = 1, nfl1
c	 delphi(ii) = delphi(ii) + delint*ttt(ii)
         if (inarea(ii).eq.1) then
	   rx(ii) = rx(ii) + delint*ttx(ii)/ttt(ii)
	   ry(ii) = ry(ii) + delint*tty(ii)/ttt(ii)
	   rz(ii) = rz(ii) + delint*ttz(ii)/ttt(ii)
c	 if (tty(ii).le.0.1) write(*,*) 'tty<0', ii,tty(ii),
c     +                                  rx(ii),ry(ii),rz(ii)
        endif
  200 continue
c
      return
      end
      subroutine bintarc
c***********************************************************************
	include 'misflin'
c
      integer   ii, indfast, zpmnull(nfl1), inarf(nfl1)
      integer   iz0,iz1,iv00,iv01,iv10,iv11
      integer   ix(nfl1), iy(nfl1), iz(nfl1), iv(nfl1),
     +          ixf(nfl1), iyf(nfl1), izf(nfl1)
      real      hxu(nfl1), hyu(nfl1), hzu(nfl1),
     +          hxo(nfl1), hyo(nfl1), hzo(nfl1)
      real      hxf(nfl1), hyf(nfl1), hzf(nfl1),
     +          hxf1(nfl1), hyf1(nfl1), hzf1(nfl1)
      real      ttx(nfl1), tty(nfl1), ttz(nfl1)
      real      ttt(nfl1)
      real      rxa(nfl1), rya(nfl1), rza(nfl1)
c
c  Integration of field lines using a 2 step 
c           lax wendroff scheme.
c   >indfast< determines grid indices for field interpolation.
c   Care is required for special symmetries!
c   Symmetries considered are 
c   (a) linesymmetry along x axis for mpause (mpause=true)
c   (b) linesymmetry along y axis for arc (mpause=false)
c   Use of symmetry requires mirror=true
c.........................................
c I. Step
c
      do 20 ii = 1, nfl1
	if (rx(ii).ge.x3d(nx-1))   inarea(ii) = 0
	if (rx(ii).le.x3d(   2))   inarea(ii) = 0
	if (ry(ii).ge.y3d(ny-1))   inarea(ii) = 0
	if (ry(ii).le.y3d(   2))   inarea(ii) = 0
	if (rz(ii).ge.z3d(nz-1))   inarea(ii) = 0
	if (mirror .and. rz(ii).le.-z3d(nz-1))  inarea(ii) = 0
	if (.not.mirror .and. rz(ii).le.z3d(2))  inarea(ii) = 0
   20 continue
      isum=0
      do 22 ii = 1, nfl1
        isum=isum+inarea(ii)
   22 continue
      if (isum.eq.0) then
        inatot=.false.
        return
      endif

      do 30 ii = 1, nfl1
	if (rz(ii) .lt. 0) then
	  zpmnull(ii) = -1
	else
	  zpmnull(ii) = 1
	end if
   30 continue

      do 50 ii = 1, nfl1
       if (inarea(ii).eq.1) then
	ix(ii) = indfast (zpmnull(ii)*rx(ii),x3d,nx)
	iy(ii) = indfast (ry(ii),y3d,ny)
        iz(ii) = indfast (zpmnull(ii)*rz(ii),z3d,nz)
        iv(ii) = ix(ii) + (iy(ii)-1)*nx
	hxo(ii) = (zpmnull(ii)*rx(ii)-x3d(ix(ii))) / dx3d(ix(ii))
        hyo(ii) = (ry(ii)-y3d(iy(ii))) / dy3d(iy(ii))
        hzo(ii) = (zpmnull(ii)*rz(ii)-z3d(iz(ii))) / dz3d(iz(ii))
       endif
   50 continue
      do 60 ii = 1, nfl1
       if (inarea(ii).eq.1) then
	 hxu(ii) = 1.0 - hxo(ii)
	 hyu(ii) = 1.0 - hyo(ii)
	 hzu(ii) = 1.0 - hzo(ii)
       endif
   60 continue
c
      do 90 ii = 1, nfl1
       if (inarea(ii).eq.1) then
	 ix0 = ix(ii)
	 ix1 = ix(ii) + 1
	 iy0 = iy(ii)
	 iy1 = iy(ii) + 1
	 iz0 = iz(ii)
	 iz1 = iz(ii) + 1
	 ttx(ii) = zpmnull(ii)*(
     +               hxu(ii) * hyu(ii) * hzu(ii) * bx(ix0,iy0,iz0)
     +             + hxu(ii) * hyu(ii) * hzo(ii) * bx(ix0,iy0,iz1)
     +             + hxu(ii) * hyo(ii) * hzu(ii) * bx(ix0,iy1,iz0)
     +             + hxu(ii) * hyo(ii) * hzo(ii) * bx(ix0,iy1,iz1)
     +             + hxo(ii) * hyu(ii) * hzu(ii) * bx(ix1,iy0,iz0)
     +             + hxo(ii) * hyu(ii) * hzo(ii) * bx(ix1,iy0,iz1)
     +             + hxo(ii) * hyo(ii) * hzu(ii) * bx(ix1,iy1,iz0)
     +             + hxo(ii) * hyo(ii) * hzo(ii) * bx(ix1,iy1,iz1) )
	 tty(ii) = ( hxu(ii) * hyu(ii) * hzu(ii) * by(ix0,iy0,iz0)
     +             + hxu(ii) * hyu(ii) * hzo(ii) * by(ix0,iy0,iz1)
     +             + hxu(ii) * hyo(ii) * hzu(ii) * by(ix0,iy1,iz0)
     +             + hxu(ii) * hyo(ii) * hzo(ii) * by(ix0,iy1,iz1)
     +             + hxo(ii) * hyu(ii) * hzu(ii) * by(ix1,iy0,iz0)
     +             + hxo(ii) * hyu(ii) * hzo(ii) * by(ix1,iy0,iz1)
     +             + hxo(ii) * hyo(ii) * hzu(ii) * by(ix1,iy1,iz0)
     +             + hxo(ii) * hyo(ii) * hzo(ii) * by(ix1,iy1,iz1) )
	 ttz(ii) = zpmnull(ii)*( 
     +             hxu(ii) * hyu(ii) * hzu(ii) * bz(ix0,iy0,iz0)
     +             + hxu(ii) * hyu(ii) * hzo(ii) * bz(ix0,iy0,iz1)
     +             + hxu(ii) * hyo(ii) * hzu(ii) * bz(ix0,iy1,iz0)
     +             + hxu(ii) * hyo(ii) * hzo(ii) * bz(ix0,iy1,iz1)
     +             + hxo(ii) * hyu(ii) * hzu(ii) * bz(ix1,iy0,iz0)
     +             + hxo(ii) * hyu(ii) * hzo(ii) * bz(ix1,iy0,iz1)
     +             + hxo(ii) * hyo(ii) * hzu(ii) * bz(ix1,iy1,iz0)
     +             + hxo(ii) * hyo(ii) * hzo(ii) * bz(ix1,iy1,iz1) )
       endif
  90  continue
c
      do 95 ii = 1, nfl1
       if (inarea(ii).eq.1) then
         ttt(ii)=sqrt(ttx(ii)*ttx(ii)+tty(ii)*tty(ii)+ttz(ii)*ttz(ii))
       endif
  95  continue
      do 100 ii = 1, nfl1
       if (inarea(ii).eq.1) then
	rxa(ii) = rx(ii) + 0.5*delint*ttx(ii)/ttt(ii)
	rya(ii) = ry(ii) + 0.5*delint*tty(ii)/ttt(ii)
	rza(ii) = rz(ii) + 0.5*delint*ttz(ii)/ttt(ii)
       endif
  100 continue
c
c   Scalar flux field:
c
      do 103 ii = 1, nfl1
  103   inarf(ii) = 1
      do 105 ii = 1, nfl1
	if (rxa(ii).ge.xfmax .or. rxa(ii).le.xfmin) inarf(ii) = 0
	if (rya(ii).ge.yfmax .or. rya(ii).le.yfmin) inarf(ii) = 0
	if (rza(ii).ge.zfmax .or. rza(ii).le.zfmin) inarf(ii) = 0
  105 continue
      do 110 ii = 1, nfl1
       if (inarf(ii).eq.1) then
         ixf(ii)  = int((rxa(ii)-xfmin)/xfsize*(nxf-1))+1
         iyf(ii)  = int((rya(ii)-yfmin)/yfsize*(nyf-1))+1
         izf(ii)  = int((rza(ii)-zfmin)/zfsize*(nzf-1))+1
         hxf(ii)  = (rxa(ii)-xf(ixf(ii)))/delxf
         hyf(ii)  = (rya(ii)-yf(iyf(ii)))/delyf
         hzf(ii)  = (rza(ii)-zf(izf(ii)))/delzf
         hxf1(ii) = 1.-hxf(ii)
         hyf1(ii) = 1.-hyf(ii)
         hzf1(ii) = 1.-hzf(ii)
       endif
  110 continue
      do 115 ii = 1, nfl1
       if (inarf(ii).eq.1) then
         ff(ixf(ii),iyf(ii),izf(ii))=
     +     ff(ixf(ii),iyf(ii),izf(ii))+hxf1(ii)*hyf1(ii)*hzf1(ii)
         ff(ixf(ii)+1,iyf(ii),izf(ii))=
     +     ff(ixf(ii)+1,iyf(ii),izf(ii))+hxf(ii)*hyf1(ii)*hzf1(ii)
         ff(ixf(ii),iyf(ii)+1,izf(ii))=
     +     ff(ixf(ii),iyf(ii)+1,izf(ii))+hxf1(ii)*hyf(ii)*hzf1(ii)
         ff(ixf(ii)+1,iyf(ii)+1,izf(ii))=
     +     ff(ixf(ii)+1,iyf(ii)+1,izf(ii))+hxf(ii)*hyf(ii)*hzf1(ii)
         ff(ixf(ii),iyf(ii),izf(ii)+1)=
     +     ff(ixf(ii),iyf(ii),izf(ii)+1)+hxf1(ii)*hyf1(ii)*hzf(ii)
         ff(ixf(ii)+1,iyf(ii)+1,izf(ii))=
     +     ff(ixf(ii)+1,iyf(ii)+1,izf(ii))+hxf(ii)*hyf(ii)*hzf1(ii)
         ff(ixf(ii),iyf(ii)+1,izf(ii)+1)=
     +     ff(ixf(ii),iyf(ii)+1,izf(ii)+1)+hxf1(ii)*hyf(ii)*hzf(ii)
         ff(ixf(ii)+1,iyf(ii),izf(ii)+1)=
     +     ff(ixf(ii)+1,iyf(ii),izf(ii)+1)+hxf(ii)*hyf1(ii)*hzf(ii)
         ff(ixf(ii)+1,iyf(ii)+1,izf(ii)+1)=
     +     ff(ixf(ii)+1,iyf(ii)+1,izf(ii)+1)+hxf(ii)*hyf(ii)*hzf(ii)
       endif
  115 continue
c
c II. Step
c
      do 120 ii = 1, nfl1
       if (inarea(ii).eq.1) then
	if (rxa(ii).ge.x3d(nx-1))   inarea(ii) = 0
	if (rxa(ii).le.x3d(   2))   inarea(ii) = 0
	if (rya(ii).ge.y3d(ny-1))   inarea(ii) = 0
	if (rya(ii).le.y3d(   2))   inarea(ii) = 0
	if (rza(ii).ge.z3d(nz-1))   inarea(ii) = 0
	if (mirror .and. rza(ii).le.-z3d(nz-1))  inarea(ii) = 0
	if (.not.mirror .and. rza(ii).le.z3d(2))  inarea(ii) = 0
       endif
  120 continue
      isum=0
      do 122 ii = 1, nfl1
        isum=isum+inarea(ii)
  122 continue
      if (isum.eq.0) then
        inatot=.false.
        return
      endif
      do 130 ii = 1, nfl1
	 if (rza(ii) .lt. 0) then
	    zpmnull(ii) = -1
	 else
	    zpmnull(ii) = 1
	 end if
  130 continue
c
      do 150 ii = 1, nfl1
       if (inarea(ii).eq.1) then
	 ix(ii) = indfast (zpmnull(ii)*rxa(ii),x3d,nx)
	 iy(ii) = indfast (rya(ii),y3d,ny)
	 iz(ii) = indfast (zpmnull(ii)*rza(ii),z3d,nz)
	 iv(ii) = ix(ii) + (iy(ii)-1)*nx
	 hxo(ii) = (zpmnull(ii)*rxa(ii)-x3d(ix(ii))) / dx3d(ix(ii))
	 hyo(ii) = (rya(ii)-y3d(iy(ii))) / dy3d(iy(ii))
	 hzo(ii) = (zpmnull(ii)*rza(ii)-z3d(iz(ii))) / dz3d(iz(ii))
       endif
 150  continue
      do 160 ii = 1, nfl1
       if (inarea(ii).eq.1) then
	 hxu(ii) = 1.0 - hxo(ii)
	 hyu(ii) = 1.0 - hyo(ii)
	 hzu(ii) = 1.0 - hzo(ii)
       endif
 160  continue
c
      do 190 ii = 1, nfl1
       if (inarea(ii).eq.1) then
	 ix0 = ix(ii)
	 ix1 = ix(ii) + 1
	 iy0 = iy(ii)
	 iy1 = iy(ii) + 1
	 iz0 = iz(ii)
	 iz1 = iz(ii) + 1
	 ttx(ii) = zpmnull(ii)*(
     +               hxu(ii) * hyu(ii) * hzu(ii) * bx(ix0,iy0,iz0)
     +             + hxu(ii) * hyu(ii) * hzo(ii) * bx(ix0,iy0,iz1)
     +             + hxu(ii) * hyo(ii) * hzu(ii) * bx(ix0,iy1,iz0)
     +             + hxu(ii) * hyo(ii) * hzo(ii) * bx(ix0,iy1,iz1)
     +             + hxo(ii) * hyu(ii) * hzu(ii) * bx(ix1,iy0,iz0)
     +             + hxo(ii) * hyu(ii) * hzo(ii) * bx(ix1,iy0,iz1)
     +             + hxo(ii) * hyo(ii) * hzu(ii) * bx(ix1,iy1,iz0)
     +             + hxo(ii) * hyo(ii) * hzo(ii) * bx(ix1,iy1,iz1) )
	 tty(ii) = ( hxu(ii) * hyu(ii) * hzu(ii) * by(ix0,iy0,iz0)
     +             + hxu(ii) * hyu(ii) * hzo(ii) * by(ix0,iy0,iz1)
     +             + hxu(ii) * hyo(ii) * hzu(ii) * by(ix0,iy1,iz0)
     +             + hxu(ii) * hyo(ii) * hzo(ii) * by(ix0,iy1,iz1)
     +             + hxo(ii) * hyu(ii) * hzu(ii) * by(ix1,iy0,iz0)
     +             + hxo(ii) * hyu(ii) * hzo(ii) * by(ix1,iy0,iz1)
     +             + hxo(ii) * hyo(ii) * hzu(ii) * by(ix1,iy1,iz0)
     +             + hxo(ii) * hyo(ii) * hzo(ii) * by(ix1,iy1,iz1) )
	 ttz(ii) = zpmnull(ii)*( 
     +             hxu(ii) * hyu(ii) * hzu(ii) * bz(ix0,iy0,iz0)
     +             + hxu(ii) * hyu(ii) * hzo(ii) * bz(ix0,iy0,iz1)
     +             + hxu(ii) * hyo(ii) * hzu(ii) * bz(ix0,iy1,iz0)
     +             + hxu(ii) * hyo(ii) * hzo(ii) * bz(ix0,iy1,iz1)
     +             + hxo(ii) * hyu(ii) * hzu(ii) * bz(ix1,iy0,iz0)
     +             + hxo(ii) * hyu(ii) * hzo(ii) * bz(ix1,iy0,iz1)
     +             + hxo(ii) * hyo(ii) * hzu(ii) * bz(ix1,iy1,iz0)
     +             + hxo(ii) * hyo(ii) * hzo(ii) * bz(ix1,iy1,iz1) )
       endif
 190  continue
c
      do 195 ii = 1, nfl1
       if (inarea(ii).eq.1) then
         ttt(ii)=sqrt(ttx(ii)*ttx(ii)+tty(ii)*tty(ii)+ttz(ii)*ttz(ii))
       endif
 195  continue
      do 200 ii = 1, nfl1
c	 delphi(ii) = delphi(ii) + delint*ttt(ii)
	 rx(ii) = rx(ii) + delint*ttx(ii)/ttt(ii)
	 ry(ii) = ry(ii) + delint*tty(ii)/ttt(ii)
	 rz(ii) = rz(ii) + delint*ttz(ii)/ttt(ii)
  200 continue
c
      return
      end
      integer function indfast (ort,gitt,max)
c***********************************************************************
      real     ort, gitt(*)
      integer  max,l, r, k
c  determines i such that  gitt(i) <= ort <= gitt(i)  
c.................................................................
      if (ort.lt.gitt(1)) then
	indfast = 1
      else
        if (ort.gt.gitt(max)) then
	  indfast = gitt(max-1)
        else
	   l = max
	   r = 1
  10      continue
	   k = (l+r) / 2
	   if (ort.le.gitt(k)) l = k - 1
	   if (ort.ge.gitt(k)) r = k + 1
	  if (r.le.l) goto 10
	  if (r-l.eq.2) then
	    indfast = l + 1
	  else
	    indfast = l
	  endif
        endif
      endif
      return
      end
      subroutine orgjout
c***********************************************************************
	include 'misflin'
c
	character*9     mjcray
	character*26     mjoutn
c.................................................................
	do 100 ix = 1+ndiffj, nx-ndiffj
        do 100 iy = 2, ny-1
        do 100 iz = 2, nz-1
	  feldx(ix,iy,iz) = dify(iy)*(bz(ix,iy+1,iz)-bz(ix,iy-1,iz))
     +                   - difz(iz)*(by(ix,iy,iz+1)-by(ix,iy,iz-1))
	  feldy(ix,iy,iz) = difz(iz)*(bx(ix,iy,iz+1)-bx(ix,iy,iz-1))
     +                   - difx(ix)*(bz(ix+1,iy,iz)-bz(ix-1,iy,iz))
	  feldz(ix,iy,iz) = difx(ix)*(by(ix+1,iy,iz)-by(ix-1,iy,iz))
     +                   - dify(iy)*(bx(ix,iy+1,iz)-bx(ix,iy-1,iz))
	  hilf(ix,iy,iz) = feldx(ix,iy,iz)*feldx(ix,iy,iz)
     +                  + feldy(ix,iy,iz)*feldy(ix,iy,iz)
     +                  + feldz(ix,iy,iz)*feldz(ix,iy,iz)
  100   continue
	do 200 ix = 1, nxj
        do 200 iy = 1, nyj
        do 200 iz = 1, nzj
	  ffj(ix,iy,iz) = sqrt(
     +         feldx(ix+ndiffj,iy+1,iz+1)*feldx(ix+ndiffj,iy+1,iz+1)
     +       + feldy(ix+ndiffj,iy+1,iz+1)*feldy(ix+ndiffj,iy+1,iz+1)
     +       + feldz(ix+ndiffj,iy+1,iz+1)*feldz(ix+ndiffj,iy+1,iz+1)
     +                        + 0.000001 )
  200   continue

         if (ieeout .eq. 1) then
           write(mjoutn,1) njout
           call asnunit(23,mjoutn,ier)
         else
           write(mjcray,2) njout
	   open (23,file=mjcray,form='unformatted')
	 endif
	 	 
         write(23) nxj,nyj,nzj,time
         write(23) xj,yj,zj
         write(23) ffj
	 
	 if (ieeout .eq. 1)   call assign('assign -R')
         close (23)
        njout=njout+1
        
    1 format('-a mjout',i3.3,' -F f77 -N ieee')
    2 format('mjcray',i3.3)
c
      return
      end
	subroutine nzznew
c***********************************************************************
	include 'misflin'
c
	integer  ix,iy,iz,newnz,nzzold
c ........................................................
c  termination of program for rho,u .lt. 0 oder v .gt. 4
c ........................................................
      nzzold = nzz
      newnz = nzz-4
      do 30 iz = nzz-3,nzz-1
      do 30 iy = 1,ny
      do 30 ix = 1,nx
       if ( abs (rho(ix,iy,iz)-rhoprof(ix) ) .ge. cval .or.
     +      abs (u(ix,iy,iz)-uprof(ix) )     .ge. cval .or.
     +      abs (sx(ix,iy,iz)-sxprof(ix) )   .ge. cval .or.
     +      abs (sy(ix,iy,iz)-syprof(ix) )   .ge. cval .or.
     +      abs (sz(ix,iy,iz)-szprof(ix) )   .ge. cval .or.
     +      abs (bx(ix,iy,iz)-bxprof(ix) )   .ge. cval .or.
     +      abs (by(ix,iy,iz)-byprof(ix) )   .ge. cval .or.
     +      abs (bz(ix,iy,iz)-bzprof(ix) )   .ge. cval )  newnz = iz
   30 continue
      nzz = min((newnz+4),nz)
      if (nzz.gt.nzzold) write(26,1) nzz,time,z(nzz),nzzold
      if (nzz.gt.nzzold) write(*,1) nzz,time,z(nzz),nzzold
      
    1 format(' new nzz = ',i4,'   at time =',f9.4,
     +       '   z(nzz) =',f9.3,'    nzzold =',i4)
c
      return
      end
