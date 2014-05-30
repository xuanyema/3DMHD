	subroutine checkem
c***********************************************************************
	include 'misflin'
c
c       Version 2
c       This subroutine checks various variables for problematic values.
c       It was designed for use in diagnosing problems with the plasma
c       outflow.  It is compiled using Version 4 of makefile (and 
c       possibly later versions of the makefile).
c       Version 2 differs from Version 1 primarily in the formatting of
c       some output statements.
c       16 December 2002
c
	integer  ix,iz,iy

c       Generate some identifying information in the standard output.
        write(6,910) 'Entering diagnostic subroutine checkem at time ',
     +               time,'istep ',istep
c
c       Begin diagnostic loops.
        do 100 iz = 1,nz
           do 100 iy = 1,ny
              do 100 ix = 1,nx
                 if (u(ix,iy,iz) .lt. 1.0e-5) 
     +              write(6,915) 'u',ix,iy,iz,u(ix,iy,iz)
c                if (abs(bx(ix,iy,iz)) .lt. 1.0e-5) 
c    +              write(6,915) 'bx',ix,iy,iz,bx(ix,iy,iz)
c                if (abs(by(ix,iy,iz)) .lt. 1.0e-5)
c    +              write(6,915) 'by',ix,iy,iz,by(ix,iy,iz)
c                if (abs(bz(ix,iy,iz)) .lt. 1.0e-5)
c    +              write(6,915) 'bz',ix,iy,iz,bz(ix,iy,iz)
c                if (abs(sx(ix,iy,iz)) .lt. 1.0e-5)
c    +              write(6,915) 'sx',ix,iy,iz,sx(ix,iy,iz)
c                if (abs(sy(ix,iy,iz)) .lt. 1.0e-5) 
c    +              write(6,915) 'sy',ix,iy,iz,sy(ix,iy,iz)
c                if (abs(sz(ix,iy,iz)) .lt. 1.0e-5)
c    +              write(6,915) 'sz',ix,iy,iz,sz(ix,iy,iz)
c                if (rho(ix,iy,iz) .lt. 1.0e-5) 
c    +              write(6,915) 'rho',ix,iy,iz,rho(ix,iy,iz)
 100    continue
	write(6,*) 'Exiting diagnostic subroutine checkem'

c

 910    format(1x,A,f8.3,5x,', ',A,i5)
 915	format(1x,'Concern: ',A,'(',i2,',',i2,',',i2,') = ',g14.5)
c
      return
      end
 
