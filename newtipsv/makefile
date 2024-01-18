# makefile for tipsv

PROGS = mpi-tipsv tipsv
FC1 = mpif90
FC2 = ifort
FFLAGS = -O #-no-pie

SRC1 = calmat.f90 solver.f90 trialf.f90 others.f90 ../common/utils.f90 fileio.f90 mpi-tipsv.f90
SRC2 = calmat.f90 solver.f90 trialf.f90 others.f90 ../common/utils.f90 fileio.f90 tipsv.f90
OBJS1 = $(SRC1:.f90=.o) rk3.o glu2.o
OBJS2 = $(SRC2:.f90=.o) rk3.o glu2.o


.PHONY: all
all: $(PROGS)

mpi-tipsv: $(OBJS1)
	$(FC1) $(FFLAGS) -o $@ $(OBJS1)

tipsv: $(OBJS2)
	$(FC2) $(FFLAGS) -o $@ $(OBJS2)

# For compiling mpi-tipsv.o, check fileio.mod and mpi-tipsv.f90.
mpi-tipsv.o: mpi-tipsv.f90 fileio.mod
	$(FC1) $(FFLAGS) -c $< -fdec -o $@

# For compiling tipsv.o, check fileio.mod in addition to tipsv.f90.
tipsv.o: fileio.mod


rk3.o: rk3.f
	$(FC2) $(FFLAGS) -c $< -o $@

glu2.o: glu2.f
	$(FC2) $(FFLAGS) -c $< -o $@



# For .o files (excluding mpi-tipsv.o), compile when .f90 files are updated.
%.o: %.f90
	$(FC2) $(FFLAGS) -c $< -o $@

# For .mod files, try to create the .o file.
# When .f90 files are newer than .o files, they will be compiled, 
#  and .mod files may or may not be updated in the process.
%.mod: %.o
	@:

.PHONY: clean
clean:
	rm -f $(OBJS1) $(OBJS2) $(PROGS) fileio.mod
