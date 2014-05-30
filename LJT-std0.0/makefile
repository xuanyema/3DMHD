# makefile for 3D MHD code
# Version 3
# This version includes entries which allow for the inclusion
# of the Tsyganenko semi-empirical magnetic field model and
# the accompanying GeoPack suite of subroutines.
# circa 15 April 2002



OBJ =	magnfl.o	first.o		orgstart.o	\
	intstart.o	satstart.o	satdat.o	\
	finterp.o	orgint.o	osmooth.o	\
	pread.o		pwrite.o	grid.o		\
	ltinitcon.o	smooth.o	meanlax.o	\
	lax.o		leap.o		fbound0.o	\
	bdcoef.o	bound.o		resbd.o		\
	termin.o	diag0.o		diag1.o		\
	diag2.o		diag3.o		diagout.o	\
	resist.o	binin.o  	binout.o	\
	nzznew.o	fndiag.o	boundxt.o	

#FC = ifort -O1 -c
#FC = pgf90 -fast -mp -c
#FC = ftn -fast -mp=nonuma -c
#FC = gfortran -O3 -fopenmp -c
FC = gfortran -O3  -c

DEF = misflin

mess3d: $(OBJ) $(OBJT)  $(OBJG)
#	ifort -no-ipo $(OBJ) $(OBJT) $(OBJG) -o mess3d
#	pgf90 $(OBJ) -fast -mp -o mhd3d
#	ftn $(OBJ) -fast -mp=nonuma -o mhd3d
#	gfortran $(OBJ) -O3 -fopenmp -o mhd3d
	gfortran $(OBJ) -O3  -o mhd3d


magnfl.o:	magnfl.f $(DEF)		; $(FC)   magnfl.f   	
first.o:	first.f $(DEF)		; $(FC)   first.f   
orgstart.o:	orgstart.f $(DEF)	; $(FC)   orgstart.f   
intstart.o:	intstart.f $(DEF)	; $(FC)   intstart.f   
satstart.o:	satstart.f $(DEF)	; $(FC)   satstart.f   
satdat.o:	satdat.f $(DEF)		; $(FC)   satdat.f   
finterp.o:	finterp.f $(DEF)	; $(FC)   finterp.f   
orgint.o:	orgint.f $(DEF)		; $(FC)   orgint.f   
osmooth.o:	osmooth.f $(DEF)	; $(FC)   osmooth.f   
pread.o:	pread.f $(DEF)		; $(FC)   pread.f   
pwrite.o:	pwrite.f $(DEF)		; $(FC)   pwrite.f   
grid.o:		grid.f $(DEF)		; $(FC)   grid.f   	
ltinitcon.o:	ltinitcon.f $(DEF)	; $(FC)   ltinitcon.f   
smooth.o:	smooth.f $(DEF)		; $(FC)   smooth.f   	
meanlax.o:	meanlax.f $(DEF)	; $(FC)   meanlax.f   
lax.o:		lax.f $(DEF)		; $(FC)   lax.f   	
leap.o:		leap.f $(DEF)		; $(FC)   leap.f   
fbound0.o:	fbound0.f $(DEF)	; $(FC)	  fbound0.f 
bdcoef.o:	bdcoef.f $(DEF)		; $(FC)   bdcoef.f   
bound.o:	bound.f $(DEF)		; $(FC)   bound.f   	
boundxt.o:	boundxt.f $(DEF)	; $(FC)   boundxt.f   	
resbd.o:	resbd.f $(DEF)		; $(FC)   resbd.f   	
termin.o:	termin.f $(DEF)		; $(FC)   termin.f  
diag0.o:	diag0.f $(DEF)		; $(FC)   diag0.f   	
diag1.o:	diag1.f $(DEF)		; $(FC)   diag1.f   	
diag2.o:	diag2.f $(DEF)		; $(FC)   diag2.f   	
diag3.o:	diag3.f $(DEF)		; $(FC)   diag3.f   	
diagout.o:	diagout.f $(DEF)	; $(FC)   diagout.f   
resist.o:	resist.f $(DEF)		; $(FC)   resist.f   	
binin.o:	binin.f $(DEF)		; $(FC)   binin.f   	
binout.o:	binout.f $(DEF)		; $(FC)   binout.f   	
nzznew.o:	nzznew.f $(DEF)		; $(FC)   nzznew.f
fndiag.o:	fndiag.f $(DEF)		; $(FC)	  fndiag.f 

cleano:
	rm -f *.o

cleant:
	rm -f $(OBJT)

cleang:
	rm -f $(OBJG)
