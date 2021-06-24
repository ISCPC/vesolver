SRCDIR = ../src

CC = ncc
CFLAGS = -O3 -fpic -I$(SRCDIR)/common -DSXAT -D_TIMELOG
LDFLAGS = -fpic -shared -lcblas -lblas_openmp -lsblas_openmp

CSRCS = $(SRCDIR)/common/LinearSolver.c $(SRCDIR)/common/Matrix.c
CSRCS += $(SRCDIR)/solvers/bicgstab2.c
CSRCS += $(SRCDIR)/solvers/cg.c
CSRCS += $(SRCDIR)/solvers/hs.c

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

COBJS = $(CSRCS:.c=.o)
FOBJS = $(FSRCS:.c=.o)
OBJS = $(COBJS) $(FOBJS)

OBJS2 = LinearSolver.o Matrix.o bicgstab2.o cg.o hs.o

libsolver.so: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS2) $(LDFLAGS)

clean:
	rm $(OBJS2) libsolver.so
