# makefile for tipsv

PROGS = mpi-tipsv tipsv
FC = mpif90
FC2 = ifort
FFLAGS = -O -no-pie
## The following may be needed for Ubuntu (set -L option in $LDFLAGS to path to intel compiler library)
#LDFLAGS = -lifport -lifcoremt -limf -lintlc -lsvml -lirng -L /opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin

SRC = calmat.f trialf.f others.f rk3.f glu2.f dcsymbdl0.f dcsymbdl.f mpi-tipsv.f formpi.f
SRC2 = calmat.f trialf.f others.f rk3.f glu2.f dcsymbdl0.f dcsymbdl.f tipsv.f
OBJS = $(SRC:.f=.o)
OBJS2 = $(SRC2:.f=.o)
.SUFFIXES: .f .o

all: $(PROGS)

mpi-tipsv: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDFLAGS)

tipsv: $(OBJS2)
	$(FC2) $(FFLAGS) -o $@ $(OBJS2)

mpi-tipsv.o: mpi-tipsv.f
	$(FC) $(FFLAGS) -c mpi-tipsv.f -fdec -o $@

.f.o:
	$(FC2) $(FFLAGS) -c $< 

clean:
	rm -f $(OBJS) $(OBJS2) $(PROGS) work
