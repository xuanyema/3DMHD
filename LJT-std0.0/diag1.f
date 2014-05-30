 	subroutine diag1
c***********************************************************************
	include 'misflin'
c
c       Version 6p
c       This version differs only slightly from Version 5.  Output
c       pertaining to additional elements of the vmax array is
c       generated.
c       Modifications introduced by Antonius Otto and Fred Hall IV
c       30 October 2002
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

c     Write output to magdmax, other diagnostic files, and
c     the standard output.
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

c     Write time and maximum velocity magnitude to magsvmax.
      write(41,23) time,vmax(4,1),vmax(4,2),vmax(4,3),vmax(4,4)

c     Write this information to the standard output, as well.
      write(6,35)
      write(6,*) 'Diagnostic output generated'
      write(6,31) 'istep = ',istep
      write(6,33) 'Time recorded: ',time
      write(6,37) 'Time to more decimal places: ',time
      write(6,31) 'ldiag = ',ldiag
      write(6,28) vmax(4,1)
      write(6,29) vmax(4,2),vmax(4,3),vmax(4,4)

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
   23 format(2f14.5,3f14.3)
   27 format(/1x,' Time: ',f8.2)
   28 format( 1x,' Maximum Velocity Magnitude: ',f14.5)
   29 format( 1x,' Location: ',/1x,'(x,y,z) = (',
     +        f8.2,',',f8.2,',',f8.2,')')
   31 format( 1x,A,i5)
   33 format( 1x,A,f8.2)
   35 format((/))
   37 format( 1x,A,f14.5)
      end if
c
      return
      end
