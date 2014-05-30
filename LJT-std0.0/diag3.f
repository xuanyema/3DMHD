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
