OBJS = main.o InitFuncs.o RunSimulation.o OutputFuncs.o

###################################################
#   Compilers    #
###################################################
COMP=/usr/bin/g++
COMP_MPI=mpic++
CFLAGS=-O3
NAME=cru3d
INCLUDES=-I. -Ixercesc-dist/include -I/mnt/galaxyTools/tools/mpich/mpich2-1.4.1p1/include
LIBS=-lxerces-c
LIB_PATHS=-Lxercesc-dist/lib
LINKER_OPTS=-Wl,-rpath=xercesc-dist/lib

##########################################
#               main program             #
##########################################
all:	$(OBJS)
		$(COMP_MPI) $(OBJS) $(CFLAGS) $(LIBS) $(LIB_PATHS) $(LINKER_OPTS) -o $(NAME)

clean:
	rm -f *.o $(NAME)

###########################################
#               subroutines               #
###########################################

main.o: main.cpp  parameters.h
	$(COMP_MPI) -c $(CFLAGS) $(INCLUDES) $<

InitFuncs.o: InitFuncs.cpp  parameters.h
	$(COMP_MPI) -c $(CFLAGS) $(INCLUDES) $<

RunSimulation.o: RunSimulation.cpp  parameters.h
	$(COMP) -c $(CFLAGS) $(INCLUDES) $<
	
OutputFuncs.o: OutputFuncs.cpp  parameters.h
	$(COMP_MPI) -c $(CFLAGS) $(INCLUDES) $<
	
