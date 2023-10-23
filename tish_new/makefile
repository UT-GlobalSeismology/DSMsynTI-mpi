# makefile for tish

PROGS = mpi-tish tish
FC1 = mpif90
FC2 = gfortran
FFLAGS = -O #-no-pie

SRC1 = calmat.f90 solver.f90 trialf.f90 others.f90 fileio.f90 formpi.f90 mpi-tish.f90
SRC2 = calmat.f90 solver.f90 trialf.f90 others.f90 fileio.f90 tish.f90
OBJS1 = $(SRC1:.f90=.o)
OBJS2 = $(SRC2:.f90=.o)


.PHONY: all
all: $(PROGS)

mpi-tish: $(OBJS1)
	$(FC1) $(FFLAGS) -o $@ $(OBJS1)

tish: $(OBJS2)
	$(FC2) $(FFLAGS) -o $@ $(OBJS2)

# For compiling mpi-tish.o, check fileio.mod and mpi-tish.f90.
mpi-tish.o: mpi-tish.f90 fileio.mod
	$(FC1) $(FFLAGS) -c mpi-tish.f90 -fdec -o $@

# For compiling tish.o, check fileio.mod in addition to tish.f90.
tish.o: fileio.mod

# For .o files (excluding mpi-tish.o), compile when .f90 files are updated.
%.o: %.f90
	$(FC2) $(FFLAGS) -c $< 

# For .mod files, try to create the .o file.
# When .f90 files are newer than .o files, they will be compiled, 
#  and .mod files may or may not be updated in the process.
%.mod: %.o
	@:

.PHONY: clean
clean:
	rm -f $(OBJS1) $(OBJS2) $(PROGS) *.mod
