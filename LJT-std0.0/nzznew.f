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
