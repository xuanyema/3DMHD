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
