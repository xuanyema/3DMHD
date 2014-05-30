	subroutine fbound0
c***********************************************************************
	include 'misflin'
c
	integer  ix,iy,iz
c .........................................................
c  set boundary of feld to 0
c .........................................................
	 do 121 iy = 1,ny
	 do 121 ix = 1,nx
	     feldx(ix,iy,1) = 0.
	     feldy(ix,iy,1) = 0.
	     feldz(ix,iy,1) = 0.
	     feldx(ix,iy,nz) = 0.
	     feldy(ix,iy,nz) = 0.
	     feldz(ix,iy,nz) = 0.
  121    continue
      if ( perio(2) )  then
	 do 122 iz = 1,nz
	   do 122 ix = 1,nx
	     feldx(ix,1,iz) = feldx(ix,ny-2,iz)
	     feldy(ix,1,iz) = feldy(ix,ny-2,iz)
	     feldz(ix,1,iz) = feldz(ix,ny-2,iz)
	     feldx(ix,ny,iz) = feldx(ix,3,iz)
	     feldy(ix,ny,iz) = feldy(ix,3,iz)
	     feldz(ix,ny,iz) = feldz(ix,3,iz)
  122    continue
      else
	 do 123 iz = 1,nz
	   do 123 ix = 1,nx
	     feldx(ix,1,iz) = 0.
	     feldy(ix,1,iz) = 0.
	     feldz(ix,1,iz) = 0.
	     feldx(ix,ny,iz) = 0.
	     feldy(ix,ny,iz) = 0.
	     feldz(ix,ny,iz) = 0.
 123	  continue
       endif
	 do 124 iz = 1,nz
	   do 124 iy = 1,ny
	     feldx(1,iy,iz) = 0.
	     feldy(1,iy,iz) = 0.
	     feldz(1,iy,iz) = 0.
	     feldx(nx,iy,iz) = 0.
	     feldy(nx,iy,iz) = 0.
	     feldz(nx,iy,iz) = 0.
  124    continue
      return
      end
