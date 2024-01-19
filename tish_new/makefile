# makefile for tish

PROGS = mpi-tish tish
FC1 = mpif90
FC2 = ifort
FFLAGS = -O #-no-pie

OBJS0 = calmat.o solver.o trialf.o others.o operation.o utils.o fileio.o 
OBJS1 = mpi-tish.o
OBJS2 = tish.o


.PHONY: all
all: $(PROGS)

mpi-tish: $(OBJS0) $(OBJS1)
	$(FC1) $(FFLAGS) -o $@ $(OBJS0) $(OBJS1)

tish: $(OBJS0) $(OBJS2)
	$(FC2) $(FFLAGS) -o $@ $(OBJS0) $(OBJS2)

# For compiling mpi-tish.o, check fileio.mod and mpi-tish.f90.
mpi-tish.o: mpi-tish.f90 fileio.mod
	$(FC1) $(FFLAGS) -c $< -fdec -o $@

# For compiling tish.o, check fileio.mod in addition to tish.f90.
tish.o: fileio.mod

trialf.o: ../common/trialf.f90
	$(FC2) $(FFLAGS) -c $< -o $@

utils.o: ../common/utils.f90
	$(FC2) $(FFLAGS) -c $< -o $@

fileio.o: ../common/fileio.f90
	$(FC2) $(FFLAGS) -c $< -o $@

# For .o files (excluding mpi-tish.o and those in ../common/), compile when .f90 files are updated.
%.o: %.f90
	$(FC2) $(FFLAGS) -c $< -o $@

# For .mod files, try to create the .o file.
# When .f90 files are newer than .o files, they will be compiled, 
#  and .mod files may or may not be updated in the process.
%.mod: %.o
	@:

.PHONY: clean
clean:
	rm -f $(OBJS0) $(OBJS1) $(OBJS2) $(PROGS) fileio.mod
