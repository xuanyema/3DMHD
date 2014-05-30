	program magnfl
c	Version 8
c	This version is based on Version 7.  It differs in that new 
c	profiles are *not* stored at the end of the run.  The boundary
c       values at xmin and xmax are not stored either.
c       Note that the writing of the binary file at the end of the 
c       run can be turned on or off through the use of the logical
c       variable writelast.
c	Modifications introduced by Fred Hall IV
c	29 March 2004
c***********************************************************************
	include 'misflin'
c
	integer          ksafe,kdia,msafe,mout
	integer          ix
	logical          writelast
c................................................................
c  Note: 1. To run program with output of current density or 
c           flux surface set parameters in magin, modify makefile to 
c           compile additional subroutines, and uncomment lines following 
c           markers  current! or flux! in: magnfl.f, orgstart.f, misflin
c                                          
c        2. To run program on PC or workstation set ieein/ieeout to false 
c           and comment lines containing asnunit or assign in: 
c           binin.f, binout.f, and satstart.f
c................................................................
c   i) bookkeeping of program, open und close, formatted output
c................................................................
        open (15,file='magin')
        open (17,file='magdia1')
        open (18,file='magdia2')
        open (19,file='magdfte')
        open (25,file='magdmax')
        open (26,file='magpar')
        open (27,file='magdfl1')
        open (28,file='magdfl2')
        open (41,file='magsvmax')
        open (43,file='magsfnorm')
c
c       files im programm:
c       14=magtapn, 15=magin,   
c       17=magdia1, 18=magdia2, 19=magdfte, 
c       25=magdmax, 26=magpar, 
c       27=magdfl1, 28=magdfl2, 
c       31=sat1,    32=sat2,    33=sat3, 
c       34=sat4,    35=sat5,    36=sat6,    39=satbin
c       41=magsvmax
c       43=magsfnorm
c
c       Note that unit 6 corresponds to the standard output.
c
	writelast = .false.
        call orgstart
        if (fini) stop
        msafe=0
        mout=0
        ldiag = ldiag + 1
        tdiag(ldiag) = time
        call diag0
        call diag1
        call fndiag
        call diag2
c
c
c
c       Uncomment the following lines for 
c       flux! computation or
c       current! density output
c       if (time .ge. 50.) then
c          if ( fluxcal )  call orgflux
c          if ( jout )  call orgjout
c       endif

c    -----------------------------------
c   ii)  start of integration (first 2 steps)
c    -----------------------------------
 400	call intstart
        if (nzz .lt. nz) call nzznew
        if ( ferror ) go to 900
        if ( istep .ge. iend) goto 900

c    -----------------------------------
c     iii) diagnostics and output
c    -----------------------------------
 500	if (mod(istep,100).eq.0) write(26,58) istep
c    -----------------------------------
c     a) diagnostics 
c    -----------------------------------
        kdia = mod(istep,idiag)
        if ( kdia .eq. 0 )  then
	   if (ldiag .ge. ndiag) then
	      call diagout
	      ldiag = 0
	   endif
	   ldiag = ldiag + 1
	   tdiag(ldiag) = time
	   call diag0
	   call diag1
	   call fndiag
	   call diag2
        endif
        if (mod(istep,10*idiag) .eq. 0) call diag3


c    -----------------------------------
c   flux integration and current density output
c    -----------------------------------
c   Uncomment the following lines for 
c   flux! computation or
c   current! density output
c       if (time .ge. 50.) then
c          if ( fluxcal .and. (tflux .ge. deltfl) )  then
c             call orgflux
c             tflux = tflux-deltfl
c          endif
c          tflux = tflux+2.*dt
c          if ( jout .and. (tjout .ge. deltj) )  then
c             call orgjout
c             tjout = tjout-deltj
c          endif
c          tjout = tjout+2.*dt
c       endif


c    -----------------------------------
c     b) satellites
c    -----------------------------------
        if ( time .ge. tsat0 .and. mitsat)  then
           if (.not. startsat) call satstart
           call satdat
        endif


c    -----------------------------------
c     c) binary output (ksafe)
c    -----------------------------------
        ksafe = mod(istep,isafe)
        if ( ksafe.eq.0 .or. newdt) then
           igrid=1
	   call lax
	   call bound
	   call termin
	   if (newdt) dt = 0.9*dt
	   newdt = .false.
	   if ( ferror .and. iferror.le.nsmooth ) call osmooth
	   if ( ferror ) go to 900
	   if ( ksafe .eq. 0 ) then
	      msafe=msafe+1
	      call binout(msafe)
	      write(26,52) istep
	   endif
        endif
        if ( istep .ge. iend) goto 900
        if ( ksafe.eq.0 ) goto 400
c
        call orgint
        if (nzz .lt. nz) call nzznew
        if ( ferror ) goto 900
c
        goto 500
 900	call diagout
        if (writelast) then
c
c	   Write the last binary file at the end of the run.
	   msafe=msafe+1
	   call binout(msafe)
	   write(26,52) istep
	endif
c

        close (15)
        close (17)
        close (18)
        close (19)
        close (25)
        close (26)
        close (27)
        close (28)
        close (41)
        close (43)
        if ( startsat)  then
c          call assign('assign -R')
           close (49)
           close (31)
           close (32)
           close (33)
           close (34)
           close (35)
           close (36)
        endif
c
c
c
   1    format(1x,l10)
  52    format(1x,'binary data written',i5)
  58    format(1x,'at integration step:',i6)
c
c       list of subroutines:
c       ------------------------
c       main magn:
c       block data vorbes0:
c             orgstart:
c             intstart:
c             orgint:
c             oglatt:
c             lespar:
c             schreib:
c             grid:
c             anfang:
c             orgaus:
c             aus:
c             glatt:
c             meanlax:
c             lax:
c             leap:
c             rdcoef:
c             rand:
c             symfx:
c             symafx:
c             symtx:
c             symty:
c             symtz:
c             abbruch:
c             diag1:
c             fndiag:
c             diag2:
c             diagaus:
c
        stop
        end
