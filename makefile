# makefile for program GranPrix

OBJ =  main.o rw_conf17.o init3.o  findforces3.o  smoothNonExtended.o updatePressure5.o  \
    rotate.o  generate3.o  linklist3.o  force.o  vel_ver3.o  pssubs6.o   
   

SRC = main.f90 rw_conf17.f90  init3.f90 findforces3.f90  smoothNonExtended.f90 updatePressure5.f90 \
     rotate.f90  generate3.f90 linklist3.f90 force.f90  vel_ver3.f90  pssubs6.f90    
	
MODULES_SRC = mycommons.f90


FFLAGS  =  -m64 -O3 -ffixed-line-length-0 -ffree-line-length-0 -Wunused-parameter -Wunderflow  -Wline-truncation -fimplicit-none -Waliasing -fbounds-check
FC = gfortran

full: parser config OBJ GranFrixrm

GranFrixrm: OBJ  
	gfortran $(FFLAGS) -OS -o GranFrixrm $(OBJ) config_template.o mycommons.o xmlparse.a $(LIBxgotoS) 

OBJ:  MODULES config_template.o xmlparse.a
	gfortran $(FFLAGS) -c $(SRC) 
	
MODULES:
	gfortran $(FFLAGS) -c $(MODULES_SRC)
	
config: config/config_template.xml
	make -C config
	
parser:
	make -C xml-fortran/src

clean: 
	rm -f *.lst $(OBJ) GranFrixrm

fullclean:	clean 
	rm -f *.mod *.a *.o *.so
	

