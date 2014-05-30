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
           write(*,*) 'Change the values of ieein and ieeout in magin' 
c          call asnunit(14,comagt,ier)
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
	  if (ieein) write(*,*) 'Change values of ieein and ieeout' 
c       call assign('assign -R')
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
