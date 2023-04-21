# makefile for tipsv
# For CentOS, you don't have to set $(LDFLAGS) and -fdec option
# For Ubuntu, you have to set -L option in $(LDFLAGS) to path to intel compiler library
PROGS = mpi-tipsv
FC = mpif90
FC2 = ifort
FFLAGS = -O -no-pie
LDFLAGS = -lifport -lifcoremt -limf -lintlc -lsvml -lirng -L /opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin

SRC = calmat.f trialf.f others.f rk3.f glu2.f dcsymbdl0.f dcsymbdl.f mpi-tipsv.f formpi.f
SRC2 = calmat.f trialf.f others.f rk3.f glu2.f dcsymbdl0.f dcsymbdl.f tipsv.f
OBJS = $(SRC:.f=.o)
OBJS2 = $(SRC2:.f=.o)
.SUFFIXES: .f .o

all:$(PROGS) tipsv

mpi-tipsv: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDFLAGS)

tipsv: $(OBJS2)
	$(FC2) $(FFLAGS) -o $@ $(OBJS2)

mpi-tipsv.o: mpi-tipsv.f
	$(FC) $(FFLAGS) -c mpi-tipsv.f -fdec -o $@

.f.o:
	$(FC2) $(FFLAGS) -c $< 

clean:
	rm -f $(OBJS) $(OBJS2) $(PROGS) tipsv work
