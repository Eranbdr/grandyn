# Naive makefile for the simple test program "tstparse"
#
# $Id: makefile,v 1.6 2006/06/09 10:30:30 arjenmarkus Exp $
#
# Coarse test for the XML-Fortran module
#
# Note:
# You may need to change the name of the Fortran compiler
# and the compile options
#
#FC	=	f90
#FC	=	lf90
#FC	=	ifort
#FOPT	=	-c
FC	=	gfortran
FOPT	=	-m64 -O3 -c

all		:	xmlparse.a xmlreader install

install		:	xmlparse.a xmlreader
	cp *.a ../../config/
	cp *.mod ../../config/
	cp xmlreader ../../config/
	
	

xmlparse.a	:	xmlparse.o read_xml_prims.o write_xml_prims.o
	ar r xmlparse.a xmlparse.o read_xml_prims.o write_xml_prims.o

xmlparse.o	:	xmlparse.f90
	$(FC) $(FOPT) xmlparse.f90

#tstparse	:	tstparse.f90 xmlparse.a #write_xml_prims.o
#	$(FC) $(FOPT) tstparse.f90
#	$(FC) -o tstparse tstparse.o xmlparse.a

xmlreader	:	xmlreader.f90 xmlparse.a
	$(FC) $(FOPT) xmlreader.f90
	$(FC) -o xmlreader xmlreader.o xmlparse.a

read_xml_prims.o	:	read_xml_prims.f90 xmlparse.f90
	$(FC) $(FOPT) read_xml_prims.f90

write_xml_prims.o	:	write_xml_prims.f90 xmlparse.f90
	$(FC) $(FOPT)  write_xml_prims.f90

clean:	
	rm -f *.o
	rm -f *.mod
	rm -f *.a
