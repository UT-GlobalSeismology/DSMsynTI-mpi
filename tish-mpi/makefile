# makefile for tish

PROGS = mpi-tish tish
FC = mpif90
FC2 = ifort
FFLAGS = -O -no-pie
## The following may be needed for Ubuntu (set -L option in $LDFLAGS to path to intel compiler library)
#LDFLAGS = -lifport -lifcoremt -limf -lintlc -lsvml -lirng -L /opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin

SRC = calmat.f trialf.f others.f dclisb.f dclisb3.f formpi.f mpi-tish.f
SRC2 = calmat.f trialf.f others.f dclisb.f dclisb3.f tish.f
OBJS = $(SRC:.f=.o)
OBJS2 = $(SRC2:.f=.o)
.SUFFIXES: .f .o

all: $(PROGS)

mpi-tish: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDFLAGS)

tish: $(OBJS2)
	$(FC2) $(FFLAGS) -o $@ $(OBJS2)

mpi-tish.o: mpi-tish.f
	$(FC) $(FFLAGS) -c mpi-tish.f -fdec -o $@

.f.o:
	$(FC2) $(FFLAGS) -c $< 

clean:
	rm -f $(OBJS) $(OBJS2) $(PROGS) work
