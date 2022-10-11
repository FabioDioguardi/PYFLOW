# fortran compiler
F90 = gfortran
# bin directory
CODE_BIN = .
#COMPFLAG = -w -O3 -msse2 -fdefault-real-8 -ffast-math 
#COMPFLAG = -w -O0 

#-------------------------------------------------------------------------------
# compiling an linking options
FFLAGS = -c $(COMPFLAG)
#FLD = -v $(COMPFLAG)
# rules
.SUFFIXES:
.SUFFIXES: .f .f90 .o

default: pyflow 

OBJS = nrtype.o nrutil.o inputdata.o\
	fminln.o \

ROUT = testt.f90 qtrap.f90 zbrent.f90\
	newt.f90 functions.f90 probfunction.f90\
	profiles.f90 twocomponent.f90\
	twolayer.f90 main.f90 

pyflow: $(OBJS) 
	$(F90) -o $(CODE_BIN)/pyflow.exe $(OBJS) $(ROUT) 

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	del *.o *.mod
