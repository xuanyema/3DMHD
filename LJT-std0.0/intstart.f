	subroutine intstart
c***********************************************************************
	include 'misflin'
c............................................................
c       integration der mhd-gleichungen:
c          1. integr. step = lax step on grid 1  
c          (igrid=1).  leap step for lax wendroff always on 
c          grid 0 for the leapfrog alternating on both grids
c............................................................
c
c    first two integrationssteps:
c...............................................
c
      write(26,61) istep
   61 format(1x,'output after initcon, integration step',i5)
c      call orgout
c      call symtx

      igrid = 1
      call lax
      call bound

      istep = istep + 1
      time = time + dt
      write(26,62) istep
   62 format(1x,'output after laxstep, integration step',i5)
c      call orgout
c      call symtx

      call termin
      if ( ferror .and. iferror.le.nsmooth ) call osmooth
      if ( ferror ) return
c
      igrid = 0
      call leap
      call bound
c
      istep = istep + 1
      time = time + dt
      write(26,63) istep
   63 format(1x,'output after leapstep, integration step',i5)
c      call orgout
c      call symtx

      call termin
      if ( ferror .and. iferror.le.nsmooth ) call osmooth
c
c      write(6,60) istep
c   60 format(1x,'integration step',i5)
c
      return
      end
