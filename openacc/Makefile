SHELL='bash'
#
# Bridges - PSC
#
# Intel Compilers are loaded by default
# You will need to specifically switch to GNU Modules
# With with `modules.sh`
#

CC = g++
PGCPP = pgc++
LIBS =


TARGETS = openacc

all:	$(TARGETS)

openacc: openacc.o common.o
	$(PGCPP) -acc -o $@ $(LIBS) openacc.o common.o

openacc.o: openacc.cpp common.h
	$(PGCPP) -c -acc -Minfo=accel  openacc.cpp

common.o: common.cpp common.h
	$(PGCPP) -c -acc common.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
