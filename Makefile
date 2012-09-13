.SUFFIXES: .cxx .o

EXPAND = Spherical
DEVICE = CPU

CXX = g++ -ggdb3 -Wall -Wextra -Wshadow -Wuninitialized -O3 -ffast-math -fopenmp -funroll-loops -fforce-addr -fPIC -I./include

LFLAGS = -D$(DEVICE) -D$(EXPAND)

OBJECT = kernel/$(DEVICE)$(EXPAND)Laplace.o kernel/CPUP2P.o

serial: test/serialrun.cxx $(OBJECT)
	$(CXX) $? $(LFLAGS) -o serial

.cxx.o:
	$(CXX) -c $? -o $@ $(LFLAGS)
