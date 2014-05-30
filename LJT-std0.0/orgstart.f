	subroutine orgstart
c***********************************************************************
	include 'misflin'
c
c       Version 9p
c       Version 9p differs from Version 8p only in the addition of the
c       line initializing the logical variable ndtprt.
c       Modified by Fred Hall IV
c       2 February 2004
c
c..........................................
c    startprozeduren:
c         -parameter read
c         -boundconditions
c         -initial state
c..........................................

      write(17,34)
      write(18,35)
      write(27,36)
      write(28,37)
      write(25,38)
c
      rewind 15
      read (15,1)  start,ieein,ieeout,relax,jout,fluxcal,mitsat
      write(26,999)
  999 format(1x,'start read')
      ndtprt = .true.
c    -----------------------------------------------------------
c     i) start of program by 
c            a: analytical state (start = true)
c            b: data from magtap (start = false)
c    -----------------------------------------------------------
      if (start) then
	 call pread
	 iferror = 0
	 istep = 0
	 time = 0.
	 call grid
	 call bdcoef
	 call initcon
	 ferror = .false.
  150    continue
	 if ( (zentr(1)) .and. dt .lt. 0.5/difx((nx+1)/2) .or.
     +        (.not.(zentr(1))) .and. dt .lt. 0.5/difx(2) ) go to 200
	 dt = .5*dt
	 go to 150
  200    continue
      else
         call binin
	 call pread
	 call bdcoef
	 call initcon
      end if
      if (nzz.gt.nz) nzz=nz
      call pwrite
      istep = 0
      
      startsat=.false.
      if (mitsat .and. tsat0 .eq. 0.) then 
       call satstart
       call satdat
      endif
c
c     Apply the boundary conditions.
      call bound
c
c uncomment the following lines for 
c   flux! computation 
c  or current! density output
c      if (fluxcal)  call initflux
c  grid for jout
c      do 910 ix = 1,nxj
c  910     xj(ix)=x(ix+ndiffj)
c      do 920 iy = 1,nyj
c  920     yj(iy)=y(iy+1)
c      do 930 iz = 1,nzj
c  930     zj(iz)=z(iz+1)
c      write(*,*) 'xj boundaries:', xj(1),xj(nxj)
c      write(*,*) 'yj boundaries:', yj(1),yj(nyj)
c      write(*,*) 'zj boundaries:', zj(1),zj(nzj)
  
c
    1 format(1x,l10)
    2 format(1x,i9)
    3 format(1x,f14.5)
c
   34 format(/1x,' time  !   masse    !  masse msp !',
     +           '  kin ener  !  mag ener  !  the energ',
     +       /1x,'-------!------------!------------!',
     +           '------------!------------!------------')
c
   35 format(/1x,' time  !v grad p/2!v dot jxb !res mal j2!',
     +                    'd/dt ekin ! d/dt eb  ! d/dt eth ',
     +       /1x,'-------!----------!----------!----------!',
     +                    '----------!----------!----------')
c
   36 format(/1x,' time  !   fluss    ! flussaend. !',
     +           ' aend ideal ! aend resist',
     +       /1x,'-------!------------!------------!',
     +           '------------!--------------')
c
   37 format(/1x,' time  !flux: mass! kin. en. ! the. en. !',
     +                    '  b. en.  !ideal.Con.! res. Con.',
     +       /1x,'-------!----------!----------!----------!',
     +                    '----------!----------!----------')
c
   38 format(/1x,' TIME  ',
     +       /1x,' MAXIMUM !   X   !    Y  !    Z  !',
     +           ' MINIMUM !   X   !    Y  !    Z   ',
     +       /1x,'---------!-------!-------!-------!',
     +           '---------!-------!-------!--------')
c
      return
      end
