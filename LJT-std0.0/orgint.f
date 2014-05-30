	subroutine orgint
c***********************************************************************
	include 'misflin'
c................................................................
c
c  1)  leapfrog integration:
c.................................
c
      if (intart .eq. 1) then
	 igrid = 1
	 call leap
	 call bound
	 istep = istep + 1
	 time = time + dt
c	 call symtx
	 call termin
	 if ( ferror .and. iferror.le.nsmooth ) call osmooth
	 if ( ferror ) return

	 igrid = 0
	 call leap
	 call bound
	 istep = istep + 1
	 time = time + dt
	 if (mod(istep,100).eq.0) write(6,60) istep
c	 call symtx
	 call termin
	 if ( ferror .and. iferror.le.nsmooth ) call osmooth
      end if
c
c
c  2)  lax-wendroff integration:
c.....................................
c
      if (intart .eq. 2) then
	 igrid = 1
	 call meanlax
	 call bound
c	 call symtx
	 call lax
	 call bound
	 istep = istep + 1
	 time = time + dt
c	 call symtx
	 call termin
	 if ( ferror .and. iferror.le.nsmooth ) call osmooth
	 if ( ferror ) return

	 igrid = 0
	 call leap
	 call bound
	 istep = istep + 1
	 time = time + dt
c	 write(6,60) istep
c 	 call symtx
	 call termin
	 if ( ferror .and. iferror.le.nsmooth ) call osmooth
      end if
c
   60 format(1x,'integration step',i5)
c
      return
      end
