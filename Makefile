####################################################################
#                                                                  #
#             Makefile for FIRE Solver Benchmarks                  #
#                                                                  #
####################################################################

CC = mpicc
CFLAGS = -Wall -g -O0 
CFLAGS += $(METIS_INC)
LIBS = -lm -lmpi -lmetis 
LIBS += $(METIS_LIB)
LIBPATH = '../usr/lib/openmpi/lib/libmpi.a'

INCLUDEPATH = '../usr/lib/openmpi/include/mpi.h'

LIBPOS=libpos.a
AR = ar
ARFLAGS = rv

SRCS = initialization.c compute_solution.c finalization.c test_functions.c util_read_files.c util_write_files.c
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

all: gccg 

%.o: %.c $(DEPS)
	$(CC) -I$(INCLUDEPATH) -c -o $@ $< $(CFLAGS)

gccg: gccg.c $(LIBPOS)
	$(CC) -o $@ $^ $(CFLAGS) -L$(LIBPATH) $(LIBS) 

$(LIBPOS) : $(OBJS)
	$(AR) $(ARFLAGS) $@ $+

clean:
	rm -rf *.o gccg $(LIBPOS)
