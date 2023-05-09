# fortran compiler
F90 = gfortran
# bin directory
CODE_BIN = .
#COMPFLAG = -w -O3 -msse2 -fdefault-real-8 -ffast-math 
#COMPFLAG = -g -Wall -Wtabs -ffpe-trap=invalid,zero,overflow,underflow,inexact,denormal #-fcheck=all# 
#-------------------------------------------------------------------------------
# compiling an linking options
FFLAGS = -c $(COMPFLAG)
#FLD = -v $(COMPFLAG)
# rules
.SUFFIXES:
.SUFFIXES: .f .f90 .o

default: pyflow_2.4

OBJS = nrtype.o nrutil.o inoutdata.o\
	fminln.o nr.o\

ROUT = init_variables.f90 read_data.f90 check_data.f90\
       write_data_summary.f90 functions.f90 testt.f90\
       qsimp.f90 zbrent.f90 rtflsp.f90 funcv.f90 newt.f90\
       cd_calculator.f90 testchi.f90 grainsize.f90\
       twolayer.f90 twocomponent.f90  profiles.f90\
       probfunction.f90 deposition.f90\
       write_results.f90 main.f90 

pyflow_2.4: $(OBJS) 
	$(F90) $(COMPFLAG) -o $(CODE_BIN)/pyflow_2.4.exe $(OBJS) $(ROUT) 

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	del *.o *.mod *~f90
