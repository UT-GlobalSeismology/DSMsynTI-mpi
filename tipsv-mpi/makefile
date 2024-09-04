# makefile for tipsv

PROGS = mpi-tipsv tipsv
FC1 = mpif90
FC2 = ifort
FFLAGS = -O #-no-pie

OBJS0 = calmat.o solver.o trialf.o others.o operation.o rk3.o glu2.o utils.o fileio.o
OBJS1 = mpi-tipsv.o
OBJS2 = tipsv.o


tipsv: $(OBJS0) $(OBJS2)
	$(FC2) $(FFLAGS) -o $@ $(OBJS0) $(OBJS2)

.PHONY: all
all: $(PROGS)

mpi-tipsv: $(OBJS0) $(OBJS1)
	$(FC1) $(FFLAGS) -o $@ $(OBJS0) $(OBJS1)

# For compiling mpi-tipsv.o, check fileio.mod and mpi-tipsv.f90.
mpi-tipsv.o: mpi-tipsv.f90 fileio.mod
	$(FC1) $(FFLAGS) -c $< -fdec -o $@

# For compiling tipsv.o, check fileio.mod in addition to tipsv.f90.
tipsv.o: fileio.mod

trialf.o: ../common/trialf.f90
	$(FC2) $(FFLAGS) -c $< -o $@

utils.o: ../common/utils.f90
	$(FC2) $(FFLAGS) -c $< -o $@

fileio.o: ../common/fileio.f90
	$(FC2) $(FFLAGS) -c $< -o $@

# For .o files (excluding mpi-tipsv.o and those in ../common/), compile when .f90 files are updated.
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
