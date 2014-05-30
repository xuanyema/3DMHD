	block data first
c***********************************************************************
	include 'misflin'
c.................................................................
      data    mass/ndiag*0./,massu/ndiag*0./,
     +        enmag/ndiag*0./,enkin/ndiag*0./, enthe/ndiag*0./,
     +        vgradp/ndiag*0./, vdjxb/ndiag*0./, resj2/ndiag*0./,
     +        ekinpu/ndiag*0./, ethpu/ndiag*0./, ebpu/ndiag*0./,
     +        frho/ndiag*0./, fekin/ndiag*0./, fethe/ndiag*0./,
     +        feb/ndiag*0./, febi/ndiag*0./, febr/ndiag*0./,
     +        fphi/ndiag*0./, fphip/ndiag*0./, fphipi/ndiag*0./,
     +        fphipr/ndiag*0./
      data    isat/40/,timesat/0.0/,
     +        xpos/npos*0./,ypos/npos*0./,zpos/npos*0./,
     +        dxpos/npos*0./,dypos/npos*0./,dzpos/npos*0./,
     +        vxpos/npos*0./,vypos/npos*0./,vzpos/npos*0./,
     +        xsat/nsat*0./,ysat/nsat*0./,zsat/nsat*0./,
     +        vxsat/nsat*0./,vysat/nsat*0./,vzsat/nsat*0./,
     +        bxs/nsat*0./,bys/nsat*0./,bzs/nsat*0./,
     +        vxs/nsat*0./,vys/nsat*0./,vzs/nsat*0./,
     +        rhosat/nsat*0./,psat/nsat*0./,bsat/nsat*0./,
     +        ptotsat/nsat*0./
      data    fini/.false./,start/.true./,etasw/0/,
     +        mitsat/.true./,ferror/.false./,newdt/.false./,
     +        istep/0/,ieein/.true./,isafe/800/,iend/8000/,intart/1/,
     +        ivisrho/1200/,ivisu/40000/,nsmooth/80/,ieeout/.true./,
     +        iferror/0/,igrid/1/,
     +        ldiag/0/,idiag/40/,tdiag/ndiag*0./,
     +        dt/.025/,time/0./,gamma/1.667/,
     +        eta/0.02/,visx/0.02/,visy/0.02/,visz/0.02/
c uncomment the following lines for 
c   flux! computation or current! density output
c      data    tflux/0.0/,deltfl/1.0/,nflux/0/,
c     +        tjout/0.0/,deltj/1.0/,njout/0/
      end
