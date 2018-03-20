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

OBJ_NADVIBS2MCTDH = constants.o \
	global.o \
	n2m_global.o \
	channels.o \
	iomod.o \
	utils.o \
	nadvibs2mctdh.o

OBJ_PLTNADVIBS = constants.o \
	global.o \
	n2m_global.o \
	channels.o \
	iomod.o \
	utils.o \
	pltmod.o \
	pltnadvibs.o

qc2vcham: $(OBJ)
	$(F90) $(F90OPTS) $(OBJ) $(LIBS) -o qc2vcham
	rm -f *.o *~ *.mod 2>/dev/null

pltlvc: $(OBJ_PLTLVC)
	$(F90) $(F90OPTS) $(OBJ_PLTLVC) $(LIBS) -o pltlvc
	rm -f *.o *~ *.mod 2>/dev/null

transw0: $(OBJ_TRANSW0)
	$(F90) $(F90OPTS) $(OBJ_TRANSW0) $(LIBS) -o transw0
	rm -f *.o *~ *.mod 2>/dev/null

nadvibs2mctdh: $(OBJ_NADVIBS2MCTDH)
	$(F90) $(F90OPTS) $(OBJ_NADVIBS2MCTDH) $(LIBS) -o nadvibs2mctdh
	rm -f *.o *~ *.mod 2>/dev/null

pltnadvibs: $(OBJ_PLTNADVIBS)
	$(F90) $(F90OPTS) $(OBJ_PLTNADVIBS) $(LIBS) -o pltnadvibs
	rm -f *.o *~ *.mod 2>/dev/null

%.o: %.f90
	$(F90) -c $(F90OPTS) $<

%.o: %.f
	$(F77) -c $(F77OPTS) $<

%.o: %.c
	$(CC) $(CCOPTS)  -c $<

clean_all:
	rm -f *.o *~ *.mod
