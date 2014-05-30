	subroutine binout(ind)
c***********************************************************************
	include 'misflin'
c
c       Version 2
c       A small amount of additional diagnostic code has been added
c       to the subroutine.
c       Modifications introduced by Fred Hall IV
c       17 July 2002
c
        integer          ind,ier
	character*9      magtapn
        character*27     comagt
c................................................................
         write(magtapn,1969) ind
         write(comagt,1970) ind
 1969    format('magtap',i3.3)
 1970    format('-a magtap',i3.3,' -F f77 -N ieee')
c  uncomment 'call asnunit' for cray	 
         if (ieeout) then 
	   write(*,*) 'Change the values of ieein and ieeout'
c          call asnunit(14,comagt,ier)
         else
	   open (14,file=magtapn,form='unformatted')
	 end if
	 write(6,*) 'Binary output now being generated'
	 write(6,28) 'istep = ',istep
	 write(6,29) 'Reported time = ',time
	 write(6,31) 'Time to more significant figures: time = ',
     +               time
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
	 if (ieeout) write(*,*) 'Change the values of ieein/ieeout' 
c       call assign('assign -R')
	 close (14)

   28 format(1x,A,i5)
   29 format(1x,A,f8.2)
   31 format(1x,A,f14.5)
	 
      return
      end
