AGENTS = ./agents
LIB_PATH = ./lib
CC = g++
LIB_OBJ = CTRNN.o LinearBody.o random.o TSearch.o
VPATH = lib

all: library bryden_oscillation.agt bryden_phaselag.agt lagged_sinusoid.agt linear_distance.agt clean
	@ echo Look in ./build/

# Agents

%.agt:
	@ echo Building $@
	@ cd $(AGENTS)/$@; $(CC) -o $@ -I../../lib *.cpp $(addprefix ../../lib/, $(LIB_OBJ))
	@ mkdir -p build
	@ mv $(AGENTS)/$@/$@ build

# Library

library: 
	@ echo Building library
	@ cd $(LIB_PATH); $(CC) -c *.cpp


clean:
	@ echo Cleaning object files.
	@ rm $(LIB_PATH)/*.o
	
