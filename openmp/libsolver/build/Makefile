DISTSDIR=../../dists

SRCDIR = ../src
DESTDIR = $(DISTSDIR)/ve/lib
INCLUDEDIR = $(DISTSDIR)/ve/include

CC = ncc
CFLAGS = -O3 -fpic -I$(SRCDIR)/common -DSXAT -D_TIMELOG
LDFLAGS = -fpic -shared -lcblas -lblas_openmp -lsblas_openmp

CSRCS = $(SRCDIR)/common/LinearSolver.c $(SRCDIR)/common/Matrix.c
CSRCS += $(SRCDIR)/solvers/elmer/bicgstab2.c
CSRCS += $(SRCDIR)/solvers/calculix/cg.c
CSRCS += $(SRCDIR)/solvers/heterosolver/hs.c

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

install: libsolver.so
	cp -i libsolver.so $(DESTDIR)
	cp -i $(SRCDIR)/common/LinearSolver.h $(SRCDIR)/common/Matrix.h $(SRCDIR)/common/PluginAPI.h $(INCLUDEDIR)
