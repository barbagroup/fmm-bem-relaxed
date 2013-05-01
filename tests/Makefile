#
# 'make'        build executable file
# 'make clean'  removes all .o and executable files
#
include Makefile.inc

# dependency directory
DEPSDIR := $(shell mkdir -p .deps; echo .deps)

# get the shell name to determine the OS
UNAME := $(shell uname)

# define the C compiler to use
CC := gcc
ifeq ($(UNAME), Linux)
CXX := g++-4.7 -std=gnu++0x
endif
ifeq ($(UNAME), Darwin)
CXX := $(shell for i in 4.7 4.6 4.5; do if g++-mp-$$i -v >/dev/null 2>&1; then echo g++-mp-$$i; exit; fi; done; echo false) -std=gnu++0x
OBJC := gcc
endif
LINK := $(CXX)

# define any compile-time flags
CFLAGS := -fopenmp -funroll-loops -O3 -W -Wall -Wextra
ifeq ($(DEBUG),1)
CFLAGS += -g -fno-inline -DDEBUG
endif
ifeq ($(PROFILE),1)
CFLAGS += -g -pg
endif
DEPCFLAGS = -MD -MF $(DEPSDIR)/$*.d -MP

# Other in-code flags
CFLAGS +=

# define any directories containing header files other than /usr/include
#   include directories like -Ipath/to/files
INCLUDES = -I. -I$(FMMW_DIR)/include -I$(FMMW_DIR)/kernel

# define any libraries to link into executable
#   To link in libraries (libXXX.so or libXXX.a) use -lXXX options
LIBS += -lm

##################
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
##################

# 'make' - default rule
EXECS += multi_level
EXECS += scaling
EXECS += correctness
EXECS += dual_correctness
EXECS += single_level
EXECS += single_level_stresslet
EXECS += multi_level_stresslet

all: $(EXECS)


multi_level: multi_level.o
	$(LINK) $(CFLAGS) $(LDFLAGS) $(LIBS) -o $@ $^

scaling: scaling.o
	$(LINK) $(CFLAGS) $(LDFLAGS) $(LIBS) -o $@ $^

correctness: correctness.o
	$(LINK) $(CFLAGS) $(LDFLAGS) $(LIBS) -o $@ $^

dual_correctness: dual_correctness.o
	$(LINK) $(CFLAGS) $(LDFLAGS) $(LIBS) -o $@ $^

single_level: single_level.o
	$(LINK) $(CFLAGS) $(LDFLAGS) $(LIBS) -o $@ $^

single_level_stresslet: single_level_stresslet.o
	$(LINK) $(CFLAGS) $(LDFLAGS) -o $@ $^

multi_level_stresslet: multi_level_stresslet.o
	$(LINK) $(CFLAGS) $(LDFLAGS) -o $@ $^


# suffix replacement rule for building .o's from .cpp's
#   $<: the name of the prereq of the rule (a .cpp file)
#   $@: the name of the target of the rule (a .o file)
%.o: %.cpp
	$(CXX) $(CFLAGS) $(DEPCFLAGS) $(INCLUDES) -c -o $@ $<

# 'make clean' - deletes all .o and temp files, exec, and dependency file
clean:
	-$(RM) *.o *~ */*~
	-$(RM) $(EXECS)
	$(RM) -r $(DEPSDIR)

DEPFILES := $(wildcard $(DEPSDIR)/*.d) $(wildcard $(DEPSDIR)/*/*.d)
ifneq ($(DEPFILES),)
-include $(DEPFILES)
endif

# define rules that do not actually generate the corresponding file
.PHONY: clean all
