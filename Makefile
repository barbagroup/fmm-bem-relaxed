.SUFFIXES: .cpp .o

CXX = g++ -std=c++0x -g -Wshadow -Wuninitialized -O3 -ffast-math -fopenmp -funroll-loops -fforce-addr -fPIC -I./include -I./kernel -I./test #-fno-inline

# Debug CXX line
#CXX = g++ -std=c++0x -g -O3 -ffast-math -fno-inline -I./include -I./test -I./kernel

LFLAGS =

main: serial

serial: test/serialrun.o
	$(CXX) $^ $(LFLAGS) -o serial

clean:
	@rm -f kernel/*.o test/*.o serial serial_new

.cpp.o:
	$(CXX) -c $? -o $@ $(LFLAGS)
