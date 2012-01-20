LIB_PATH = ./lib
CC = g++
LIB_OBJ = CTRNN.o LinearBody.o random.o TSearch.o
VPATH = lib

all: library bryden_oscillation.sim bryden_phaselag.sim lagged_sinusoid.sim linear_distance.sim clean
	@ echo Look in ./build/

# Simulations

%.sim:
	@ echo Building $@
	@ cd ./simulations/$@; $(CC) -o $@ -I../../lib *.cpp $(addprefix ../../lib/, $(LIB_OBJ))
	@ mkdir -p build
	@ mv ./simulations/$@/$@ build

# Library

library: 
	@ echo Building library
	@ cd $(LIB_PATH); $(CC) -c *.cpp


clean:
	@ echo Cleaning object files.
	@ rm $(LIB_PATH)/*.o
	
