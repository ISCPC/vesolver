DISTSDIR=../../dists

INSTALL = /usr/bin/install

SRCDIR = ../src
DESTDIR = $(DISTSDIR)/ve/lib
INCLUDEDIR = $(DISTSDIR)/ve/include

CC = ncc
CFLAGS = -O3 -fpic -I$(SRCDIR)/common -DSXAT -D_TIMELOG
LDFLAGS = -fpic -shared -lcblas -lblas_openmp -lsblas_openmp

CSRCS = $(SRCDIR)/common/LinearSolver.c
CSRCS += $(SRCDIR)/common/Matrix.c $(SRCDIR)/common/Matrix_SXAT.c
CSRCS += $(SRCDIR)/solvers/elmer/bicgstab2.c
CSRCS += $(SRCDIR)/solvers/elmer/elmer_cg.c
CSRCS += $(SRCDIR)/solvers/calculix/cg.c
CSRCS += $(SRCDIR)/solvers/heterosolver/hs.c
#CSRCS += $(SRCDIR)/solvers/lis/lis.c

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

COBJS = $(CSRCS:.c=.o)
FOBJS = $(FSRCS:.c=.o)
OBJS = $(COBJS) $(FOBJS)

OBJS2 = LinearSolver.o Matrix.o Matrix_SXAT.o bicgstab2.o elmer_cg.o cg.o hs.o
#OBJS2 += lis.o

libsolver.so: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS2) $(LDFLAGS)

clean:
	rm $(OBJS2) libsolver.so

install: libsolver.so
	$(INSTALL) -D libsolver.so $(DESTDIR)/libsolver.so
	$(INSTALL) -D $(SRCDIR)/common/LinearSolver.h $(INCLUDEDIR)/LinearSolver.h
	$(INSTALL) -D $(SRCDIR)/common/Matrix.h $(INCLUDEDIR)/Matrix.h
	$(INSTALL) -D $(SRCDIR)/common/PluginAPI.h $(INCLUDEDIR)/PluginAPI.h
