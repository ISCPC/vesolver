DESTDIR=../../dists

INSTALL = /usr/bin/install

INCLUDEDIR=../include
SRCDIR=../src

CC = gcc
CFLAGS = -std=c99 -O3 -g -fpic -DSXAT -I../include -I/opt/nec/ve/veos/include -c ../src/vesolver.c

LDFLAGS = -L/opt/nec/ve/veos/lib64 -Wl,-rpath=/opt/nec/ve/veos/lib64 -lveo

CSRCS = $(SRCDIR)/vesolver.c

.c.o :
	$(CC) $(CFLAGS) -c $<

OBJS = $(CSRCS:.c=.o)
OBJS2 = vesolver.o
LIB = libvesolver_api.so

all: $(LIB)

$(LIB): $(OBJS)
	$(CC) -shared -o $@ $(OBJS2) $(LDFLAGS)

clean:
	-rm *.o $(LIB)

install:
	$(INSTALL) -D $(LIB) $(DESTDIR)/lib/$(LIB)
	$(INSTALL) -D $(INCLUDEDIR)/vesolver.h $(DESTDIR)/include/vesolver.h
