.SUFFIXES: .cxx .o

EXPAND = Spherical
DEVICE = CPU

#CXX = g++ -g -Wall -Wextra -Wshadow -Wuninitialized -O3 -ffast-math -fopenmp -funroll-loops -fforce-addr -fPIC -I./include -I./include_new
CXX = g++ -g -Wall -Wextra -I./include_new -DTREECODE

LFLAGS = -D$(DEVICE) -D$(EXPAND)

OBJECT = kernel/$(DEVICE)$(EXPAND)Laplace.o kernel/CPUP2P.o

main: serial_new

serial: test/serialrun.o $(OBJECT)
	$(CXX) $^ $(LFLAGS) -o serial

serial_new: test/serial_new.o 
	$(CXX) $^ $(LFLAGS) -o serial_new

clean:
	@rm -f kernel/*.o test/*.o serial serial_new

.cxx.o:
	$(CXX) -c $? -o $@ $(LFLAGS)
