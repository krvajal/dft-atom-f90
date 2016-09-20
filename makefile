

FC = gfortran
SRC_FILES = double.f90 ode.f90  quad.f90 dft.f90  io_utils.f90 roots.f90 sample.f90
OBJ_FILES= $(patsubst %.f90,%.o,$(SRC_FILES)) 

FCFLAGS = -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace

PROGRAM = dft
all: $(PROGRAM)
	
$(PROGRAM): $(SRC_FILES)	
	$(FC) $(FCFLAGS) -o $@ $^


 %.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<

clean: 
	rm  *.o *.mod