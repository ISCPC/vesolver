DESTDIR=../../dists

INSTALL = /usr/bin/install 

SRCDIR = ../src
INCDIR = $(DESTDIR)/ve/include
LIBDIR = $(DESTDIR)/ve/lib

CC = ncc
CFLAGS = -O3 -fpic -I$(SRCDIR) -I$(INCDIR) -DSXAT -D_TIMELOG
LDFLAGS = -fpic -shared -L$(LIBDIR) -lsolver
LDFLAGS += -lcblas -lblas_openmp -lsblas_openmp
LDFLAGS += -lheterosolver_openmp -lblas_sequential

CSRCS = $(SRCDIR)/vesolver.c

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

OBJS = $(CSRCS:.c=.o)

OBJS2 = vesolver.o

libvesolver.so: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS2) $(LDFLAGS)

clean:
	rm $(OBJS2) libsolver.so

install:
	$(INSTALL) -D libvesolver.so $(LIBDIR)/libvesolver.so
