#-----------------------------------------------------------------------
# Compiler flags
#-----------------------------------------------------------------------

#
# gfortran
#
F90	= gfortran
F77	= gfortran
CC	= gcc
F90OPTS = -cpp -g -ffixed-line-length-none -ffree-line-length-none -fopenmp -O3 -fbacktrace
CCOPTS  = -g -O0

# External libraries
LIBS= -lblas -llapack

OBJ = constants.o \
	global.o \
	channels.o \
	iomod.o \
	utils.o \
	parsemod.o \
	ioqc.o \
	qc2vcham.o

OBJ_PLTLVC = constants.o \
	global.o \
	channels.o \
	iomod.o \
	pltmod.o \
	pltlvc.o

OBJ_TRANSW0 = constants.o \
	global.o \
	transmod.o \
	channels.o \
	iomod.o \
	utils.o \
	ioqc.o \
	transw0.o

qc2vcham: $(OBJ)
	$(F90) $(F90OPTS) $(OBJ) $(LIBS) -o qc2vcham
	rm -f *.o *~ *.mod 2>/dev/null

pltlvc: $(OBJ_PLTLVC)
	$(F90) $(F90OPTS) $(OBJ_PLTLVC) $(LIBS) -o pltlvc
	rm -f *.o *~ *.mod 2>/dev/null

transw0: $(OBJ_TRANSW0)
	$(F90) $(F90OPTS) $(OBJ_TRANSW0) $(LIBS) -o transw0
	rm -f *.o *~ *.mod 2>/dev/null

%.o: %.f90
	$(F90) -c $(F90OPTS) $<

%.o: %.f
	$(F77) -c $(F77OPTS) $<

%.o: %.c
	$(CC) $(CCOPTS)  -c $<

clean_all:
	rm -f *.o *~ *.mod
