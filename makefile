

FC = gfortran
SRC_FILES = double.f90 ode.f90   solvers.f90 quad.f90 roots.f90 lsda.f90 dft.f90  io_utils.f90  main.f90

UTILS_DIR = $(HOME)/fortran-utils/src/

-O -Wall -fcheck=all -g -fbacktrace

FCFLAGS = -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -I$(UTILS_DIR)   -ffpe-trap=zero,underflow,overflow
FLDFLAGS =  -llapack -L$(UTILS_DIR) -lfortran_utils

PROGRAM = dft-atom.x
all: $(PROGRAM)
	
$(PROGRAM): $(SRC_FILES)	
	$(FC) $(FCFLAGS)  $? -o dftatom.x $(FLDFLAGS)
	
 %.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<

clean: 
	rm  *.o *.mod