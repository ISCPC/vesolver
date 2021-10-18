#DESTDIR=/opt/local
DESTDIR=$(HOME)/local
#DESTDIR=$(PWD)/dists

LIBSOLVER = libsolver.so
LIBVESOLVER = libvesolver.so
LIBVESOLVER_API = libvesolver_api.so

all: $(LIBSOLVER) $(LIBVESOLVER) $(LIBSOLVER_API)

# for libsolver
$(LIBSOLVER):
	cd libsolver/build; \
	$(MAKE)

$(LIBSOLVER)-install: $(LIBSOLVER)
	cd libsolver/build; \
	$(MAKE) DESTDIR=$(DESTDIR) install

$(LIBSOLVER)-clean:
	cd libsolver/build; \
	$(MAKE) clean

# for libvesolver
$(LIBVESOLVER):
	cd libvesolver/build; \
	$(MAKE)

$(LIBVESOLVER)-install: $(LIBVESOLVER)
	cd libvesolver/build; \
	$(MAKE) DESTDIR=$(DESTDIR) install

$(LIBVESOLVER)-clean:
	cd libvesolver/build; \
	$(MAKE) clean


# for libvesolver_api
$(LIBVESOLVER_API):
	cd libvesolver_api/build; \
	$(MAKE)

$(LIBVESOLVER_API)-install: $(LIBVESOLVER_API)
	cd libvesolver_api/build; \
	$(MAKE) DESTDIR=$(DESTDIR) install

$(LIBVESOLVER_API)-clean:
	cd libvesolver_api/build; \
	$(MAKE) clean


install: $(LIBSOLVER)-install $(LIBVESOLVER)-install $(LIBVESOLVER_API)-install

clean: $(LIBSOLVER)-clean $(LIBVESOLVER)-clean $(LIBVESOLVER_API)-clean
