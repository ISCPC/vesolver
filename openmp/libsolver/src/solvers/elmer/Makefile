WITH_MKL=true
WITH_MKL_SBLAS=true
WITH_ONEAPI=false

CXXFLAGS = -g
CFLAGS = -g

BUILD_ENV = gcc

include ../standalone/Make.$(BUILD_ENV)

#######
SOLVERDIR = .

CMDS = bicgstab2 elmer_cg
SOLVER_SRC = bicgstab2.c elmer_cg.c
SOLVER_OBJ = bicgstab2.o elmer_cg.o

all: $(CMDS)

$(SOLVER_OBJ): $(SOLVER_SRC)
	$(CC) $(CFLAGS) -c $^

bicgstab2: $(COMMON_OBJS) bicgstab2.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

elmer_cg: $(COMMON_OBJS) elmer_cg.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm *.o $(CMDS)
