# makefile for tish

PROGS = mpi-tish tish
FC1 = mpif90
FC2 = gfortran
FFLAGS = -O #-no-pie

SRC1 = calmat.f90 trialf.f90 others.f90 formpi.f90 mpi-tish.f90
SRC2 = calmat.f90 trialf.f90 others.f90 tish.f90
SRC77 = dclisb.f dclisb3.f
OBJS1 = $(SRC1:.f90=.o)
OBJS2 = $(SRC2:.f90=.o)
OBJS77 = $(SRC77:.f=.o)
.SUFFIXES: .f .f90 .o

all: $(PROGS)

mpi-tish: $(OBJS1) $(OBJS77)
	$(FC1) $(FFLAGS) -o $@ $(OBJS1) $(OBJS77)

tish: $(OBJS2) $(OBJS77)
	$(FC2) $(FFLAGS) -o $@ $(OBJS2) $(OBJS77)

mpi-tish.o: mpi-tish.f90
	$(FC1) $(FFLAGS) -c mpi-tish.f90 -fdec -o $@

.f.o:
	$(FC2) $(FFLAGS) -c $< 

.f90.o:
	$(FC2) $(FFLAGS) -c $< 

clean:
	rm -f $(OBJS1) $(OBJS2) $(OBJS77) $(PROGS) *.mod
