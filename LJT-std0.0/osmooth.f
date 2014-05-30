	subroutine osmooth
c***********************************************************************
	include 'misflin'
c
c.................................................................
      iferror = iferror + 1
      write(26,1) iferror
c      write(6,1) iferror
      ferror = .false.
      call smooth
      call bound
      call termin
c
    1 format(1x,'anzahl der glaettungen: ',i3)
      return
      end
