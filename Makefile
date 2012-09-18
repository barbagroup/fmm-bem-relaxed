.SUFFIXES: .cpp .o

EXPAND = Spherical
DEVICE = CPU

CXX = g++ -ggdb3 -Wall -Wextra -Wshadow -Wuninitialized -O3 -ffast-math -fopenmp -funroll-loops -fforce-addr -fPIC -I./include

LFLAGS = -D$(DEVICE) -D$(EXPAND)

OBJECT = kernel/$(DEVICE)$(EXPAND)Laplace.o kernel/CPUP2P.o

serial: test/serialrun.o $(OBJECT)
	$(CXX) $^ $(LFLAGS) -o serial

clean:
	@rm -f kernel/*.o test/*.o serial

.cpp.o:
	$(CXX) -c $? -o $@ $(LFLAGS)
