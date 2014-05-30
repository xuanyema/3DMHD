	subroutine termin
c***********************************************************************
c
c       Version 8
c       This version differs from Version 7 in that it has been
c       formatted so that some of the logic is easier to follow upon
c       casual inspection.  Some errors related to the profiles in rho 
c       and u were also corrected.
c       The variable vcrit was replaced with v2crit to indicate that
c       a critical value of the speed squared is being used.
c       A value of v2crit more comparable to that used in the
c       earliest versions of this subroutine has been used.  That used
c       in the earlier versions incorporating local smoothing was
c       quite high.
c       Again, Version 7 is based on Version 6.  It restores the form
c       of the condition used in Version 3 and earlier forms of the
c       subroutine.  Instead of crying foul if either rho or u drop
c       below a certain fraction of their initial values, the code
c       checks if they drop below certain constant values (given by
c       rhocrit and ucrit).
c       26 February 2004.
c
	include 'misflin'
c
	integer  ix,iy,iz,ksum,tstabx(nx),tstaby(ny),tstabz(nz),
     +           ixphase,iyphase,izphase,ixgrid,iygrid,izgrid,
     +           indx(nx),indy(ny),indz(nz),is,iymax(nx),izmax(nx),
     +           nox,noy,noz
        integer  kx
	real     absv2,vphase,vpgrid,dgrid,vp(nx),rhocrit,pcrit,v2crit,
     +           ucrit,entcrit,rhoc2,gaminv,ppp1,ppp2,rhoav,uav
c ........................................................
c  termination of program for rho,u .lt. 0 oder v .gt. 4
c  date of changing: 080403:
c  correction to factor of rhocrit and ucrit
c  smoothing also around disturbed point
c ........................................................
	ppp1 =  0.25
        ppp2 =  (1.-ppp1)/3.
        gaminv = 1./gamma
        rhocrit=0.09
        rhoc2=25
        pcrit=0.14
	v2crit=400.0
	entcrit=4.0
c       v2crit=3200.
c       v2crit=9600.
        ucrit=(pcrit/2.)**gaminv

        do 20 ix = 1,nx
           tstabx(ix) = 0
   20   continue
c
        do 22 iy = 1,ny
           tstaby(iy) = 0
   22   continue
c
        do 24 iz = 1,nz
           tstabz(iz) = 0
   24   continue
c
c
c
        do 30 iz = 1,nz
        do 30 iy = 1,ny
        do 30 ix = 1,nx
c	02.10.02:
           rhocor(ix,iy,iz)=.false.
           ucor(ix,iy,iz)=.false.
           rhocor(ix,iy,iz)=.false.
           if ( (rho(ix,iy,iz) .le. rhocrit) .or.
     +          (rho(ix,iy,iz) .ge. rhoc2) .or.
     +          (u(ix,iy,iz) .le. ucrit) .or.
     +          (u(ix,iy,iz)/rho(ix,iy,iz) .ge. entcrit) .or.
     +          ( sx(ix,iy,iz)*sx(ix,iy,iz) 
     +          + sy(ix,iy,iz)*sy(ix,iy,iz)
     +          + sz(ix,iy,iz)*sz(ix,iy,iz) )/
     +           (rho(ix,iy,iz)*rho(ix,iy,iz))
     +          .ge. v2crit )  then
              write(*,*) 'correction place: ix, iy, iz:  ',ix, iy, iz
              tstabx(ix) = 1
              tstaby(iy) = 1
              tstabz(iz) = 1
c
c             Additional diagnostic code
	      write(26,63)
	      write(26,63)
	      write(26,47) ix,iy,iz,time,istep
	      write(26,48) x(ix),y(iy),z(iz)
	      write(26,49) rho(ix,iy,iz)
c              write(26,800) (rho(kx,iy-1,iz-1),kx=ix-1,ix+1)
c              write(26,801) (rho(kx,iy,iz-1),kx=ix-1,ix+1)
c              write(26,801) (rho(kx,iy+1,iz-1),kx=ix-1,ix+1)
c              write(26,802) (rho(kx,iy-1,iz),kx=ix-1,ix+1)
              write(26,802) (rho(kx,iy,iz),kx=ix-1,ix+1)
c              write(26,802) (rho(kx,iy+1,iz),kx=ix-1,ix+1)
c              write(26,803) (rho(kx,iy-1,iz+1),kx=ix-1,ix+1)
c              write(26,803) (rho(kx,iy,iz+1),kx=ix-1,ix+1)
c              write(26,803) (rho(kx,iy+1,iz+1),kx=ix-1,ix+1)

	      write(26,51) u(ix,iy,iz),2.*u(ix,iy,iz)**gamma
c              write(26,810) (2.*u(kx,iy-1,iz-1)**gamma,kx=ix-1,ix+1)
c              write(26,801) (2.*u(kx,iy,iz-1)**gamma,kx=ix-1,ix+1)
c              write(26,801) (2.*u(kx,iy+1,iz-1)**gamma,kx=ix-1,ix+1)
c              write(26,802) (2.*u(kx,iy-1,iz)**gamma,kx=ix-1,ix+1)
              write(26,802) (2.*u(kx,iy,iz)**gamma,kx=ix-1,ix+1)
c              write(26,802) (2.*u(kx,iy+1,iz)**gamma,kx=ix-1,ix+1)
c              write(26,803) (2.*u(kx,iy-1,iz+1)**gamma,kx=ix-1,ix+1)
c              write(26,803) (2.*u(kx,iy,iz+1)**gamma,kx=ix-1,ix+1)
c              write(26,803) (2.*u(kx,iy+1,iz+1)**gamma,kx=ix-1,ix+1)

c              write(26,820) (2.*u(kx,iy-1,iz-1)**gamma/
c     +                 rho(kx,iy-1,iz-1),kx=ix-1,ix+1)
c              write(26,801) (2.*u(kx,iy,iz-1)**gamma/
c     +                 rho(kx,iy,iz-1),kx=ix-1,ix+1)
c              write(26,801) (2.*u(kx,iy+1,iz-1)**gamma/
c     +                 rho(kx,iy+1,iz-1),kx=ix-1,ix+1)
c              write(26,802) (2.*u(kx,iy-1,iz)**gamma/
c     +                 rho(kx,iy-1,iz),kx=ix-1,ix+1)
              write(26,802) (2.*u(kx,iy,iz)**gamma/
     +                 rho(kx,iy,iz),kx=ix-1,ix+1)
c              write(26,802) (2.*u(kx,iy+1,iz)**gamma/
c     +                 rho(kx,iy+1,iz),kx=ix-1,ix+1)
c              write(26,803) (2.*u(kx,iy-1,iz+1)**gamma/
c     +                 rho(kx,iy-1,iz+1),kx=ix-1,ix+1)
c              write(26,803) (2.*u(kx,iy,iz+1)**gamma/
c     +                 rho(kx,iy,iz+1),kx=ix-1,ix+1)
c              write(26,803) (2.*u(kx,iy+1,iz+1)**gamma/
c     +                 rho(kx,iy+1,iz+1),kx=ix-1,ix+1)
	      write(26,59) ( sx(ix,iy,iz)*sx(ix,iy,iz) 
     +                     + sy(ix,iy,iz)*sy(ix,iy,iz)
     +                     + sz(ix,iy,iz)*sz(ix,iy,iz) )
     +                     /( rho(ix,iy,iz)*rho(ix,iy,iz) )
	      write(26,54) 'x',sx(ix,iy,iz)
	      write(26,54) 'y',sy(ix,iy,iz)
	      write(26,54) 'z',sz(ix,iy,iz)
	      write(26,55) 'x',sx(ix,iy,iz)/rho(ix,iy,iz)
	      write(26,55) 'y',sy(ix,iy,iz)/rho(ix,iy,iz)
	      write(26,55) 'z',sz(ix,iy,iz)/rho(ix,iy,iz)
	      write(26,57) 'x',bx(ix,iy,iz)
	      write(26,57) 'y',by(ix,iy,iz)
	      write(26,57) 'z',bz(ix,iy,iz)
	      write(26,63)
	      write(26,65) 'rhocrit',rhocrit
	      write(26,65) 'ucrit',ucrit
	      write(26,65) 'v2crit',v2crit
	      write(26,67) 'rhoprof',ix,rhoprof(ix)
	      write(26,67) 'uprof',ix,uprof(ix)
	      write(26,63)
	      write(26,63)
c
           endif
   30   continue
 800	format('Rho:', 3f7.3)
 801	format('    ', 3f7.3)
 802	format('     ', 3f7.3)
 803	format('      ', 3f7.3)
 810	format('  P:', 3f7.3)
 820	format('  T:', 3f7.3)

c
        ksum = 0
        do 50 ix = 1,nx
           ksum = ksum + tstabx(ix)
   50   continue
c
        if (ksum .ge. 1) then
           ferror = .true.
           write(26,38) ksum,istep
           write(*,38) ksum,istep
c
c
           is = 0
           do 60 ix = 1,nx
              is = is + tstabx(ix)
              if (tstabx(ix).eq.1) indx(is) = ix
   60      continue
           nox = is
c
           is = 0
           do 70 iy = 1,ny
              is = is + tstaby(iy)
              if (tstaby(iy).eq.1) indy(is) = iy
   70      continue
           noy = is
c
           is = 0
           do 80 iz = 1,nz
              is = is + tstabz(iz)
              if (tstabz(iz).eq.1) then
                 indz(is) = iz
                 write(*,*) 'step 2a, is, indz=', is,indz(is)
              endif
   80      continue
           noz = is
c
c
c
       	   do 100 iz = 1,noz
       	   do 100 iy = 1,noy
           do 100 ix = 1,nox
c	   02.10.02:
              absv2 = ( sx(indx(ix),indy(iy),indz(iz))
     +                 *sx(indx(ix),indy(iy),indz(iz))
     +                 +sy(indx(ix),indy(iy),indz(iz))
     +                 *sy(indx(ix),indy(iy),indz(iz))
     +                 +sz(indx(ix),indy(iy),indz(iz))
     +                 *sz(indx(ix),indy(iy),indz(iz)))
     +                 /rho(indx(ix),indy(iy),indz(iz))
     +                 /rho(indx(ix),indy(iy),indz(iz))
              if ( absv2 .ge. v2crit) then
     	         write(26,41) indx(ix),indy(iy),indz(iz),absv2,time
                 if ((indx(ix) .ne. 1) .AND.(indx(ix) .ne. nx) .AND.
     +               (indx(iy) .ne. 1) .AND.(indx(iy) .ne. ny) .AND.
     +               (indx(iz) .ne. 1) .AND.(indx(iz) .ne. nz)) then
                    sxyzcor(indx(ix)-1,indy(iy)-1,indz(iz)-1)=.true.     
                    sxyzcor(indx(ix)-1,indy(iy),indz(iz)-1)=.true.
                    sxyzcor(indx(ix)-1,indy(iy)+1,indz(iz)-1)=.true.
                    sxyzcor(indx(ix),indy(iy)-1,indz(iz)-1)=.true.     
                    sxyzcor(indx(ix),indy(iy),indz(iz)-1)=.true.
                    sxyzcor(indx(ix),indy(iy)+1,indz(iz)-1)=.true.
                    sxyzcor(indx(ix)+1,indy(iy)-1,indz(iz)-1)=.true.     
                    sxyzcor(indx(ix)+1,indy(iy),indz(iz)-1)=.true. 
                    sxyzcor(indx(ix)+1,indy(iy)+1,indz(iz)-1)=.true. 
                    sxyzcor(indx(ix)-1,indy(iy)-1,indz(iz))=.true.     
                    sxyzcor(indx(ix)-1,indy(iy),indz(iz))=.true.
                    sxyzcor(indx(ix)-1,indy(iy)+1,indz(iz))=.true.
                    sxyzcor(indx(ix),indy(iy)-1,indz(iz))=.true.     
                    sxyzcor(indx(ix),indy(iy),indz(iz))=.true.
                    sxyzcor(indx(ix),indy(iy)+1,indz(iz))=.true.
                    sxyzcor(indx(ix)+1,indy(iy)-1,indz(iz))=.true.     
                    sxyzcor(indx(ix)+1,indy(iy),indz(iz))=.true. 
                    sxyzcor(indx(ix)+1,indy(iy)+1,indz(iz))=.true.
                    sxyzcor(indx(ix)-1,indy(iy)-1,indz(iz)+1)=.true.
                    sxyzcor(indx(ix)-1,indy(iy),indz(iz)+1)=.true.
                    sxyzcor(indx(ix)-1,indy(iy)+1,indz(iz)+1)=.true.
                    sxyzcor(indx(ix),indy(iy)-1,indz(iz)+1)=.true.
                    sxyzcor(indx(ix),indy(iy),indz(iz)+1)=.true.
                    sxyzcor(indx(ix),indy(iy)+1,indz(iz)+1)=.true.
                    sxyzcor(indx(ix)+1,indy(iy)-1,indz(iz)+1)=.true.
                    sxyzcor(indx(ix)+1,indy(iy),indz(iz)+1)=.true. 
                    sxyzcor(indx(ix)+1,indy(iy)+1,indz(iz)+1)=.true.
	            write(*,45) indx(ix),indy(iy),indz(iz)
c     +                          sxyzcor(indx(ix),indy(iy),indz(iz))
                 endif
              endif
c
	      
	      if ( (rho(indx(ix),indy(iy),indz(iz))
     +             .le. rhocrit)
     +             .or. (rho(indx(ix),indy(iy),indz(iz))
     +             .ge. rhoc2) ) then
	         write(26,39) indx(ix),indy(iy),indz(iz),
     +                  rho(indx(ix),indy(iy),indz(iz)),time
c	         rho(indx(ix),indy(iy),indz(iz)) = rhoprof(indz(ix))
                 if ((indx(ix) .ne. 1) .AND. (indx(ix) .ne. nx) .AND.
     +               (indy(iy) .ne. 1) .AND. (indy(iy) .ne. ny) .AND.
     +               (indz(iz) .ne. 1) .AND. (indz(iz) .ne. nz)) then
                    rhoav = ppp1*rho(indx(ix),indy(iy),indz(iz))
     +                  + ppp2*(meanpx(indx(ix))
     +                             *rho(indx(ix)+1,indy(iy),indz(iz))
     +                        + meanmx(indx(ix))
     +                             *rho(indx(ix)-1,indy(iy),indz(iz))
     +                        + meanpy(indy(iy))
     +                             *rho(indx(ix),indy(iy)+1,indz(iz))
     +                        + meanmy(indy(iy))
     +                             *rho(indx(ix),indy(iy)-1,indz(iz)) 
     +                        + meanpz(indz(iz))
     +                             *rho(indx(ix),indy(iy),indz(iz)+1)
     +                        + meanmz(indz(iz))
     +                             *rho(indx(ix),indy(iy),indz(iz)-1) )
                    uav = ppp1*u(indx(ix),indy(iy),indz(iz))
     +                  + ppp2*(meanpx(indx(ix))
     +                             *u(indx(ix)+1,indy(iy),indz(iz))
     +                        + meanmx(indx(ix))
     +                             *u(indx(ix)-1,indy(iy),indz(iz))
     +                        + meanpy(indy(iy))
     +                             *u(indx(ix),indy(iy)+1,indz(iz))
     +                        + meanmy(indy(iy))
     +                             *u(indx(ix),indy(iy)-1,indz(iz)) 
     +                        + meanpz(indz(iz))
     +                             *u(indx(ix),indy(iy),indz(iz)+1)
     +                        + meanmz(indz(iz))
     +                             *u(indx(ix),indy(iy),indz(iz)-1) )
                    rho(indx(ix),indy(iy),indz(iz)) = rhoav
                    u(indx(ix),indy(iy),indz(iz)) = uav
c     +                  *rhoprof(indz(ix))
	         write(26,39) indx(ix),indy(iy),indz(iz),
     +                  rho(indx(ix),indy(iy),indz(iz)),time
                    rhocor(indx(ix)-1,indy(iy)-1,indz(iz)-1)=.true.     
                    rhocor(indx(ix)-1,indy(iy),indz(iz)-1)=.true.
                    rhocor(indx(ix)-1,indy(iy)+1,indz(iz)-1)=.true.
                    rhocor(indx(ix),indy(iy)-1,indz(iz)-1)=.true.     
                    rhocor(indx(ix),indy(iy),indz(iz)-1)=.true.
                    rhocor(indx(ix),indy(iy)+1,indz(iz)-1)=.true.
                    rhocor(indx(ix)+1,indy(iy)-1,indz(iz)-1)=.true.     
                    rhocor(indx(ix)+1,indy(iy),indz(iz)-1)=.true. 
                    rhocor(indx(ix)+1,indy(iy)+1,indz(iz)-1)=.true. 
                    rhocor(indx(ix)-1,indy(iy)-1,indz(iz))=.true.     
                    rhocor(indx(ix)-1,indy(iy),indz(iz))=.true.
                    rhocor(indx(ix)-1,indy(iy)+1,indz(iz))=.true.
                    rhocor(indx(ix),indy(iy)-1,indz(iz))=.true.     
                    rhocor(indx(ix),indy(iy),indz(iz))=.true.
                    rhocor(indx(ix),indy(iy)+1,indz(iz))=.true.
                    rhocor(indx(ix)+1,indy(iy)-1,indz(iz))=.true.     
                    rhocor(indx(ix)+1,indy(iy),indz(iz))=.true. 
                    rhocor(indx(ix)+1,indy(iy)+1,indz(iz))=.true.
                    rhocor(indx(ix)-1,indy(iy)-1,indz(iz)+1)=.true.     
                    rhocor(indx(ix)-1,indy(iy),indz(iz)+1)=.true.
                    rhocor(indx(ix)-1,indy(iy)+1,indz(iz)+1)=.true.
                    rhocor(indx(ix),indy(iy)-1,indz(iz)+1)=.true.     
                    rhocor(indx(ix),indy(iy),indz(iz)+1)=.true.
                    rhocor(indx(ix),indy(iy)+1,indz(iz)+1)=.true.	      
                    rhocor(indx(ix)+1,indy(iy)-1,indz(iz)+1)=.true.     
                    rhocor(indx(ix)+1,indy(iy),indz(iz)+1)=.true.
                    rhocor(indx(ix)+1,indy(iy)+1,indz(iz)+1)=.true.
                 endif
              endif
c
              if ( u(indx(ix),indy(iy),indz(iz))
     +            .le. ucrit ) then
                 if ((indx(ix) .ne. 1) .AND. (indx(ix) .ne. nx) .AND.
     +               (indy(iy) .ne. 1) .AND. (indy(iy) .ne. ny) .AND.
     +               (indz(iz) .ne. 1) .AND. (indz(iz) .ne. nz)) then
	         write(26,391) indx(ix),indy(iy),indz(iz),
     +                  u(indx(ix),indy(iy),indz(iz)),time
                    rhoav = ppp1*rho(indx(ix),indy(iy),indz(iz))
     +                  + ppp2*(meanpx(indx(ix))
     +                             *rho(indx(ix)+1,indy(iy),indz(iz))
     +                        + meanmx(indx(ix))
     +                             *rho(indx(ix)-1,indy(iy),indz(iz))
     +                        + meanpy(indy(iy))
     +                             *rho(indx(ix),indy(iy)+1,indz(iz))
     +                        + meanmy(indy(iy))
     +                             *rho(indx(ix),indy(iy)-1,indz(iz)) 
     +                        + meanpz(indz(iz))
     +                             *rho(indx(ix),indy(iy),indz(iz)+1)
     +                        + meanmz(indz(iz))
     +                             *rho(indx(ix),indy(iy),indz(iz)-1) )
                    uav = ppp1*u(indx(ix),indy(iy),indz(iz))
     +                  + ppp2*(meanpx(indx(ix))
     +                             *u(indx(ix)+1,indy(iy),indz(iz))
     +                        + meanmx(indx(ix))
     +                             *u(indx(ix)-1,indy(iy),indz(iz))
     +                        + meanpy(indy(iy))
     +                             *u(indx(ix),indy(iy)+1,indz(iz))
     +                        + meanmy(indy(iy))
     +                             *u(indx(ix),indy(iy)-1,indz(iz)) 
     +                        + meanpz(indz(iz))
     +                             *u(indx(ix),indy(iy),indz(iz)+1)
     +                        + meanmz(indz(iz))
     +                             *u(indx(ix),indy(iy),indz(iz)-1) )
                    rho(indx(ix),indy(iy),indz(iz)) = rhoav
                    u(indx(ix),indy(iy),indz(iz)) = uav
c                    u(indx(ix),indy(iy),indz(iz)) = uprof(indz(ix))
	         write(26,391) indx(ix),indy(iy),indz(iz),
     +                  u(indx(ix),indy(iy),indz(iz)),time
                    ucor(indx(ix)-1,indy(iy)-1,indz(iz)-1)=.true.     
                    ucor(indx(ix)-1,indy(iy),indz(iz)-1)=.true.
                    ucor(indx(ix)-1,indy(iy)+1,indz(iz)-1)=.true.
                    ucor(indx(ix),indy(iy)-1,indz(iz)-1)=.true.     
                    ucor(indx(ix),indy(iy),indz(iz)-1)=.true.
                    ucor(indx(ix),indy(iy)+1,indz(iz)-1)=.true.
                    ucor(indx(ix)+1,indy(iy)-1,indz(iz)-1)=.true.     
                    ucor(indx(ix)+1,indy(iy),indz(iz)-1)=.true. 
                    ucor(indx(ix)+1,indy(iy)+1,indz(iz)-1)=.true. 
                    ucor(indx(ix)-1,indy(iy)-1,indz(iz))=.true.     
                    ucor(indx(ix)-1,indy(iy),indz(iz))=.true.
                    ucor(indx(ix)-1,indy(iy)+1,indz(iz))=.true.
                    ucor(indx(ix),indy(iy)-1,indz(iz))=.true.     
                    ucor(indx(ix),indy(iy),indz(iz))=.true.
                    ucor(indx(ix),indy(iy)+1,indz(iz))=.true.
                    ucor(indx(ix)+1,indy(iy)-1,indz(iz))=.true.     
                    ucor(indx(ix)+1,indy(iy),indz(iz))=.true. 
                    ucor(indx(ix)+1,indy(iy)+1,indz(iz))=.true.
                    ucor(indx(ix)-1,indy(iy)-1,indz(iz)+1)=.true.     
                    ucor(indx(ix)-1,indy(iy),indz(iz)+1)=.true.
                    ucor(indx(ix)-1,indy(iy)+1,indz(iz)+1)=.true.
                    ucor(indx(ix),indy(iy)-1,indz(iz)+1)=.true.     
                    ucor(indx(ix),indy(iy),indz(iz)+1)=.true.
                    ucor(indx(ix),indy(iy)+1,indz(iz)+1)=.true.	      
                    ucor(indx(ix)+1,indy(iy)-1,indz(iz)+1)=.true.     
                    ucor(indx(ix)+1,indy(iy),indz(iz)+1)=.true.
                    ucor(indx(ix)+1,indy(iy)+1,indz(iz)+1)=.true.
                 endif
              endif

              if ( u(indx(ix),indy(iy),indz(iz))/
     +               rho(indx(ix),indy(iy),indz(iz))
     +            .ge. entcrit ) then
                 if ((indx(ix) .ne. 1) .AND. (indx(ix) .ne. nx) .AND.
     +               (indy(iy) .ne. 1) .AND. (indy(iy) .ne. ny) .AND.
     +               (indz(iz) .ne. 1) .AND. (indz(iz) .ne. nz)) then
	         write(26,39) indx(ix),indy(iy),indz(iz),
     +                  rho(indx(ix),indy(iy),indz(iz)),time
	         write(26,391) indx(ix),indy(iy),indz(iz),
     +                  u(indx(ix),indy(iy),indz(iz)),time
	         write(26,392) indx(ix),indy(iy),indz(iz),
     +                  u(indx(ix),indy(iy),indz(iz))/
     +                     rho(indx(ix),indy(iy),indz(iz)),time
                    rhoav = ppp1*rho(indx(ix),indy(iy),indz(iz))
     +                  + ppp2*(meanpx(indx(ix))
     +                             *rho(indx(ix)+1,indy(iy),indz(iz))
     +                        + meanmx(indx(ix))
     +                             *rho(indx(ix)-1,indy(iy),indz(iz))
     +                        + meanpy(indy(iy))
     +                             *rho(indx(ix),indy(iy)+1,indz(iz))
     +                        + meanmy(indy(iy))
     +                             *rho(indx(ix),indy(iy)-1,indz(iz)) 
     +                        + meanpz(indz(iz))
     +                             *rho(indx(ix),indy(iy),indz(iz)+1)
     +                        + meanmz(indz(iz))
     +                             *rho(indx(ix),indy(iy),indz(iz)-1) )
                    uav = ppp1*u(indx(ix),indy(iy),indz(iz))
     +                  + ppp2*(meanpx(indx(ix))
     +                             *u(indx(ix)+1,indy(iy),indz(iz))
     +                        + meanmx(indx(ix))
     +                             *u(indx(ix)-1,indy(iy),indz(iz))
     +                        + meanpy(indy(iy))
     +                             *u(indx(ix),indy(iy)+1,indz(iz))
     +                        + meanmy(indy(iy))
     +                             *u(indx(ix),indy(iy)-1,indz(iz)) 
     +                        + meanpz(indz(iz))
     +                             *u(indx(ix),indy(iy),indz(iz)+1)
     +                        + meanmz(indz(iz))
     +                             *u(indx(ix),indy(iy),indz(iz)-1) )
                    rho(indx(ix),indy(iy),indz(iz)) = rhoav
                    u(indx(ix),indy(iy),indz(iz)) = uav
c                    u(indx(ix),indy(iy),indz(iz)) = uprof(indz(ix))
	         write(26,39) indx(ix),indy(iy),indz(iz),
     +                  rho(indx(ix),indy(iy),indz(iz)),time
	         write(26,391) indx(ix),indy(iy),indz(iz),
     +                  u(indx(ix),indy(iy),indz(iz)),time
	         write(26,392) indx(ix),indy(iy),indz(iz),
     +                  u(indx(ix),indy(iy),indz(iz))/
     +                     rho(indx(ix),indy(iy),indz(iz)),time
                    ucor(indx(ix)-1,indy(iy)-1,indz(iz)-1)=.true.     
                    ucor(indx(ix)-1,indy(iy),indz(iz)-1)=.true.
                    ucor(indx(ix)-1,indy(iy)+1,indz(iz)-1)=.true.
                    ucor(indx(ix),indy(iy)-1,indz(iz)-1)=.true.     
                    ucor(indx(ix),indy(iy),indz(iz)-1)=.true.
                    ucor(indx(ix),indy(iy)+1,indz(iz)-1)=.true.
                    ucor(indx(ix)+1,indy(iy)-1,indz(iz)-1)=.true.     
                    ucor(indx(ix)+1,indy(iy),indz(iz)-1)=.true. 
                    ucor(indx(ix)+1,indy(iy)+1,indz(iz)-1)=.true. 
                    ucor(indx(ix)-1,indy(iy)-1,indz(iz))=.true.     
                    ucor(indx(ix)-1,indy(iy),indz(iz))=.true.
                    ucor(indx(ix)-1,indy(iy)+1,indz(iz))=.true.
                    ucor(indx(ix),indy(iy)-1,indz(iz))=.true.     
                    ucor(indx(ix),indy(iy),indz(iz))=.true.
                    ucor(indx(ix),indy(iy)+1,indz(iz))=.true.
                    ucor(indx(ix)+1,indy(iy)-1,indz(iz))=.true.     
                    ucor(indx(ix)+1,indy(iy),indz(iz))=.true. 
                    ucor(indx(ix)+1,indy(iy)+1,indz(iz))=.true.
                    ucor(indx(ix)-1,indy(iy)-1,indz(iz)+1)=.true.     
                    ucor(indx(ix)-1,indy(iy),indz(iz)+1)=.true.
                    ucor(indx(ix)-1,indy(iy)+1,indz(iz)+1)=.true.
                    ucor(indx(ix),indy(iy)-1,indz(iz)+1)=.true.     
                    ucor(indx(ix),indy(iy),indz(iz)+1)=.true.
                    ucor(indx(ix),indy(iy)+1,indz(iz)+1)=.true.	      
                    ucor(indx(ix)+1,indy(iy)-1,indz(iz)+1)=.true.     
                    ucor(indx(ix)+1,indy(iy),indz(iz)+1)=.true.
                    ucor(indx(ix)+1,indy(iy)+1,indz(iz)+1)=.true.
                 endif
              endif

  100      continue
c
c
c
           do 300 iz = 1,nz
           do 300 iy = 1,ny
           do 300 ix = 1,nx
              help(ix,iy,iz) = sqrt ( ( bx(ix,iy,iz)*bx(ix,iy,iz)
     +         + by(ix,iy,iz)*by(ix,iy,iz) + bz(ix,iy,iz)*bz(ix,iy,iz)
     +         + gamma*u(ix,iy,iz)**gamma ) / rho(ix,iy,iz) )
              hilf(ix,iy,iz) = max(difx(ix),dify(iy),difz(iz))
  300      continue
c
           vphase = 0.
           vpgrid = 0.
           do 350 ix = 1,nx
              vp(ix) = help(ix,1,1)
              iymax(ix) = 2
              izmax(ix) = 2
  350      continue
c
           do 360 iz = 1,nz
           do 360 iy = 2,ny 
           do 360 ix = 1,nx
              if (help(ix,iy,iz).gt.vp(ix)) then
                 vp(ix) = help(ix,iy,iz)
                 iymax(ix) = iy
                 izmax(ix) = iz
              endif
  360      continue
c
           do 370 ix = 1,nx 
              if (vp(ix).gt.vphase) then
                 vphase = vp(ix)
                 ixphase=ix
                 iyphase=iymax(ix)
                 izphase=izmax(ix) 
              endif    
  370      continue
c
           do 380 iz = 1,nz
           do 380 iy = 1,ny
           do 380 ix = 1,nx
              help(ix,iy,iz) = help(ix,iy,iz)*hilf(ix,iy,iz)
  380      continue
c
           do 390 ix = 1,nx
              vp(ix) = help(ix,1,1)
              iymax(ix) = 2
              izmax(ix) = 2
  390      continue
c
           do 400 iz = 1,nz
           do 400 iy = 2,ny
           do 400 ix = 1,nx
              if (help(ix,iy,iz).gt.vp(ix)) then
                 vp(ix) = help(ix,iy,iz)
                 iymax(ix) = iy
                 izmax(ix) = iz
              endif
  400      continue
c
           do 410 ix = 1,nx
              if (vp(ix).gt.vpgrid) then
                 vpgrid = vp(ix)
                 ixgrid=ix
                 iygrid=iymax(ix)
                 izgrid=izmax(ix)
              endif
  410      continue
c
           vpgrid=2.*vpgrid
           write(26,42) ixphase,iyphase,izphase,vphase,
     +                  ixgrid,iygrid,izgrid,vpgrid,dt
           if (dt .gt. 0.9/vpgrid) then
              newdt = .true.
              write(26,43) time, 0.9*dt
           endif
c
c
        endif
c
c
c
   38 format(/,'No of errors is greater than ',i3,'  at istep:',i5,/)
   39 format(' rho(',i3,',',i3,',',i3,') =',f12.7,'   time =',f9.4)
 391  format(' u(',i3,',',i3,',',i3,') =',f12.7,'   time =',f9.4)
 392  format(' u/rho(',i3,',',i3,',',i3,') =',f12.7,'   time =',f9.4)
   40 format(' u(',i3,',',i3,',',i3,')   =',f12.5,
     +       ' ucrit*uprof(',i3,')   =',f12.5,
     +       '   time =',f9.4)
   41 format(' abs v**2(',i3,',',i3,',',i3,') =',f12.4,
     +       '   time =',f9.4)
   42 format(' vphase(',i3,',',i3,',',i3,') =',f8.3,
     +       '   vpgrid(',i3,',',i3,',',i3,') =',f7.2,
     +       ' present dt =',f6.4)
   43 format(' time =',f10.5,'   newdt =', f7.4)
   45 format('termin-indx(ix),indx(iy),indx(iz),sxyzcor: ',
     +       3(i3))
   47 format(/' Potential problem, located at indices (',i3,', ',
     +          i3,', ',i3,'), at time ',f9.4,', istep ',i5)
   48 format(' (x,y,z): (',f9.4,', ',f9.4,', ',f9.4,')')
   49 format(' Density rho: ',f9.4)
   51 format(' Pressure analogue u and p: ',f9.4)
   53 format(' Speed (earlier calculation): ',f9.4)
   54 format(' s',A,': ',f9.4)
   55 format(' v',A,': ',f9.4)
   57 format(' b',A,': ',f9.4)
   59 format(' Speed (squared): ',f9.4)
   63 format((/))
   65 format(A,' = ',f9.4)
   67 format(A,'(',i3,') = ',f9.4)
c
      return
      end
