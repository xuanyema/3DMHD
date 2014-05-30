	subroutine resist
c***********************************************************************
c
c       Version 8
c       This version is based on Version 6.  It differs in the inclusion
c       of a new resistivity profile.  The current-dependent resistivity
c       Model 1 used in earlier versions is replaced by a new profile
c       defined in Version 37 of initcon.f.  This new profile
c       corresponds to etasw = 2.
c       20 August 2004
c       
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
	eta1 = 0.0
c       eta1 = 0.02*eta0
        do 10 iz = 1,nzz
        do 10 iy = 1,ny
        do 10 ix = 1,nx
           res(ix,iy,iz) = eta0*prof(ix,iy,iz) + eta1
   10   continue
      end if

      if (etasw .ge. 2 .and. eta.ne. 0.0) then
        jcrit1=jcrit0
c       jcrit1=jcrit0*delbz
        jcrit2=jcrit1**2
        eta0 = 0.1*eta*(exp(-time/10.0)-exp(-time/2.0))
c       eta0 = 0.0
        eta1 = 0.002
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
c    model 2
        if (etasw .ge. 3) then
          do 70 iz = 2,nzz-1
          do 70 iy = 2,ny1
	  ixanf = 2 + mod(iz+iy+igrid,2)
	  do 70 ix = ixanf,nx1,2
            if (hilf(ix,iy,iz) .gt. jcrit2)  
     +          res(ix,iy,iz) = eta*(hilf(ix,iy,iz)-jcrit2)
   70     continue
        endif
c

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
