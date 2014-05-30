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
