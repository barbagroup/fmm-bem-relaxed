.SUFFIXES: .cpp .o

CXX = g++ -g -Wshadow -Wuninitialized -O3 -ffast-math -fopenmp -funroll-loops -fforce-addr -fPIC -I./include -I./kernel #-Wall -Wextra

LFLAGS =

main: serial

serial: test/serialrun.o
	$(CXX) $^ $(LFLAGS) -o serial

clean:
	@rm -f kernel/*.o test/*.o serial serial_new

.cpp.o:
	$(CXX) -c $? -o $@ $(LFLAGS)
