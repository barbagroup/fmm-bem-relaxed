.SUFFIXES: .cpp .o

EXPAND = Spherical
DEVICE = CPU

CXX = g++ -g -Wall -Wextra -Wshadow -Wuninitialized -O3 -ffast-math -fopenmp -funroll-loops -fforce-addr -fPIC -I./include
#CXX = g++ -g -Wall -Wextra -I./include_new 

LFLAGS = -D$(DEVICE) -D$(EXPAND)

main: serial

serial: test/serialrun.o 
	$(CXX) $^ $(LFLAGS) -o serial

clean:
	@rm -f kernel/*.o test/*.o serial serial_new

.cpp.o:
	$(CXX) -c $? -o $@ $(LFLAGS)
