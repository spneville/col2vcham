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
	col2vcham.o

OBJ_PLTLVC = constants.o \
	global.o \
	channels.o \
	iomod.o \
	pltmod.o \
	pltlvc.o

col2vcham: $(OBJ)
	$(F90) $(F90OPTS) $(OBJ) $(LIBS) -o col2vcham
	rm -f *.o *~ *.mod 2>/dev/null

pltlvc: $(OBJ_PLTLVC)
	$(F90) $(F90OPTS) $(OBJ_PLTLVC) $(LIBS) -o pltlvc
	rm -f *.o *~ *.mod 2>/dev/null

%.o: %.f90
	$(F90) -c $(F90OPTS) $<

%.o: %.f
	$(F77) -c $(F77OPTS) $<

%.o: %.c
	$(CC) $(CCOPTS)  -c $<

clean_all:
	rm -f *.o *~ *.mod
