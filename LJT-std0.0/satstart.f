	subroutine satstart
c***********************************************************************
	include 'misflin'
c
      integer    i,j,k,fnum,ier
      character*8      satname
c.................................................................
c    initializing of satparameters
c.................................................................
      startsat = .true.
      timesat = time
      open (31,file='satpar')
      
      do 50 i = 1,npos
       do 50 j = 1,nrow  
         k = j + (i-1)*nrow
         xsat(k) = xpos(i) + (j-1)*dxpos(i)
         ysat(k) = ypos(i) + (j-1)*dypos(i)
         zsat(k) = zpos(i) + (j-1)*dzpos(i)
         vxsat(k) = vxpos(i)
         vysat(k) = vypos(i)
         vzsat(k) = vzpos(i)
   50 continue

      write(31,21)  
      write(26,211)  
      write(31,22) nsat,timesat,tsatout
      write(26,22) nsat,timesat,tsatout
      write(31,23) nsim,bsim,vsim,psim,lsim,tsim
      write(31,24) rhobd(1),pbd(1),vxbd(1),vybd(1),vzbd(1),
     +             bxbd(1),bybd(1),bzbd(1),
     +             rhobd(2),pbd(2),vxbd(2),vybd(2),vzbd(2),
     +             bxbd(2),bybd(2),bzbd(2)
      write(31,25) 
      write(26,25) 
      do 60 k = 1,nsat
          write(31,26) k,xsat(k),ysat(k),zsat(k),
     +                     vxsat(k),vysat(k),vzsat(k)
          write(26,26) k,xsat(k),ysat(k),zsat(k),
     +                     vxsat(k),vysat(k),vzsat(k)
   60 continue
      write(31,27) 
      write(26,27) 
   
c  uncomment for cray
      if (ieein) then 
	write(*,*) 'Change the values of ieein/ieeout'
c       call asnunit(49,'-a satbin -F f77 -N ieee',ier)
      else
        open (49,file='satbin',form='unformatted')
      end if
      write(49) nsat,timesat
      write(49) nsim,bsim,vsim,psim,lsim,tsim
      write(49) rhobd,pbd,vxbd,vybd,vzbd,bxbd,bybd,bzbd
      write(49) xsat,ysat,zsat,vxsat,vysat,vzsat

      do 70 i = 1,npos
      do 70 j = 1,nrow  
         k = j + (i-1)*nrow
         vxsat(k) = vxpos(i) + vxbd(1)
         vysat(k) = vypos(i) + vybd(1)
         vzsat(k) = vzpos(i) + vzbd(1)
   70 continue

   21 format(1x,'Satellite and program parameters:',
     +       /1x,'---------------------------------')
  211 format(1x,'Satellite parameters:',
     +       1x,'---------------------')
   22 format(/1x,'number of satellites:',i3,3x,
     +           'start time:',f7.2,3x,'time resolution:',f7.4)
   23 format(/1x,'Normalisation for : ',
     +       /1x,' No density ! magn. field!  velocity  !  pressure  !',
     +           'length units!    time     ',
     +       /1x,'  cm**(-3)  !     nT     !    km/s    !   nPascal  !',
     +           '     km     !      s      ',
     +         /,f11.2,f13.2,f13.2,f13.4,f13.2,f13.4)
   24 format(/1x,'Initial asymptotic values at xmin (1. row) ',
     +           'and xmax (2. row) for: ',
     +       /1x,'    rho      p        vx       vy       vz    ',
     +           '   bx       by       bz    ',
     +       /,8f9.3/,8f9.3)
   25 format(/1x,'Initial satellite locations and velocity:',
     +       /1x,'sat  !    x    !    y    !    z    !',
     +                 '    vx   !    vy   !    vz   ')
   26 format(1x,i3,'   ',f7.2,'   ',f7.2,'   ',f7.2,
     +             '   ',f7.2,'   ',f7.2,'   ',f7.2)
   27 format(1x,'Satellite velocity is relative to restframe of the',
     +            ' magnetosphere (at xmin)!!!',
     +       /1x,'Plasma velocity is recorded in the satellite frame!')

      return
      end
