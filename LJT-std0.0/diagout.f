	subroutine diagout
c       Version 6p
c       Version 6p is based on Version 4 of this subroutine, differing
c       only in the generation of some additional diagnostic 
c       information.
c       30 October 2002
c***********************************************************************
	include 'misflin'
c
        real rrho,rekin,rethe,reb,rebi,rebr
	integer  i
c ....................................................
c     output magdia1(17) and magdia2(18)
c     as well as a few other files
c ....................................................
      do 1100 i = 1,ldiag
            rrho = 0.0
            rekin = 0.0
            rethe = 0.0
            reb = 0.0
            rebi = 0.0
            rebr = 0.0
          if (abs(fphi(i)).gt. 0.01) then 
            rrho = frho(i)/fphi(i)
            rekin = fekin(i)/fphi(i)
            rethe = fethe(i)/fphi(i)
            reb = feb(i)/fphi(i)
            rebi = febi(i)/fphi(i)
            rebr = febr(i)/fphi(i)
          end if

	  write(17,1) tdiag(i),mass(i),massu(i),
     +                         enkin(i),enmag(i),enthe(i)
	  write(43,5) tdiag(i),fnorm(i),fnmax(i),ixmm(i),iymm(i),izmm(i)
	  write(6,7)  tdiag(i)
	  write(6,8)  fnorm(i)
	  write(18,2) tdiag(i),vgradp(i),vdjxb(i),resj2(i),
     +                ekinpu(i),ebpu(i),ethpu(i)
	  write(27,3) tdiag(i),fphi(i),fphip(i),fphipi(i),fphipr(i)
	  write(28,2) tdiag(i),frho(i),fekin(i),fethe(i),
     +                feb(i),febi(i),febr(i)
	  write(28,4) rrho,rekin,rethe,
     +                reb,rebi,rebr
 1100 continue
      do 1200 i = 1,ndiag
	  mass(i)  = 0.
	  massu(i)  = 0.
	  enmag(i) = 0.
 	  enkin(i) = 0.
	  enthe(i) = 0.
	  fnorm(i) = 0.
	  vgradp(i) = 0.
	  vdjxb(i) = 0.
	  resj2(i) = 0.
	  ekinpu(i)  = 0.
	  ethpu(i)   = 0.
	  ebpu(i)    = 0.
	  tdiag(i)   = 0.
          fphi(i)   = 0.
          fphip(i)  = 0.
          fphipi(i) = 0.
          fphipr(i) = 0.
          frho(i)   = 0.
          fekin(i)  = 0.
          fethe(i)  = 0.
          feb(i)    = 0.
          febi(i)   = 0.
          febr(i)   = 0.
 1200 continue
c
    1 format(1x,f7.3,' ',f12.2,' ',f12.2,' ',f12.3,' ',f12.2,' ',f12.2)
    2 format(1x,f7.3,' ',f10.3,' ',f10.3,' ',f10.3,' ',f10.3,
     +               ' ',f10.3,' ',f10.3 )
    3 format(1x,f7.3,' ',f12.3,' ',f12.3,' ',f12.3,' ',f12.3)
    4 format(5x,    '  ',f10.4,' ',f10.4,' ',f10.4,' ',f10.4,
     +               ' ',f10.4,' ',f10.4 )
    5 format(1x,f7.3,' ',f14.5,' ',f14.5,3i5)
    7 format(/1x,' Time: ',f8.2)
    8 format( 1x,' Force Normal: ',g8.2)

c
      return
      end
