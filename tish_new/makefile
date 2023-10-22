# makefile for tish

PROGS = mpi-tish tish
FC1 = mpif90
FC2 = gfortran
FFLAGS = -O #-no-pie

SRC1 = calmat.f90 trialf.f90 others.f90 solver.f90 formpi.f90 mpi-tish.f90
SRC2 = calmat.f90 trialf.f90 others.f90 solver.f90 tish.f90
OBJS1 = $(SRC1:.f90=.o)
OBJS2 = $(SRC2:.f90=.o)
.SUFFIXES: .f90 .o

all: $(PROGS)

mpi-tish: $(OBJS1)
	$(FC1) $(FFLAGS) -o $@ $(OBJS1)

tish: $(OBJS2)
	$(FC2) $(FFLAGS) -o $@ $(OBJS2)

mpi-tish.o: mpi-tish.f90
	$(FC1) $(FFLAGS) -c mpi-tish.f90 -fdec -o $@

.f90.o:
	$(FC2) $(FFLAGS) -c $< 

clean:
	rm -f $(OBJS1) $(OBJS2) $(PROGS) *.mod
