#############################################################################
#     Makefile for building Femocs executable from given main funciton
#############################################################################

# Flags to compile and link FEMOCS

## release mode optimization flags
FEMOCS_OPT=-O3 -w 

## debug mode optimization flags
FEMOCS_DOPT=-g -Og -Wall -Wpedantic -Wno-unused-local-typedefs -ggdb3

## paths to library
FEMOCS_PATH=lib/femocs
MD_PATH=md
MLIP_PATH=lib/mlip-2-master
DEALII_PATH=lib/deal.II
XC_PATH=lib/libxc-5.0.0

## paths to headers
HEADPATH= -I$(FEMOCS_PATH)/include -I$(FEMOCS_PATH)/lib -I$(DEALII_PATH)/include -I$(FEMOCS_PATH)/GETELEC/modules -I$(MD_PATH)/ -I$(XC_PATH)/include -I./ -std=c++14 
FEMOCS_EXT_HEAD=

## paths to libraries
LIBPATH=-L$(FEMOCS_PATH)/lib -L$(FEMOCS_PATH)/GETELEC/lib -L$(DEALII_PATH)/lib -L$(MD_PATH)/lib -L$(MLIP_PATH)/lib -L$(XC_PATH)/lib
FEMOCS_EXT_LIB=

## release mode libraries
LIB=-lfemocs -lmd -ltet -ldeal_II -lgetelec -lslatec -lxc -fopenmp -ltbb -llapack -lz -lm -lstdc++ -lgfortran -l_mlip_interface

## debug mode libraries
FEMOCS_DLIB=-lfemocs -ltet -ldeal_II -lgetelec -lslatec -lxc -fopenmp -ltbb -llapack -lz -lm -lstdc++ -lgfortran 

## flag indicating whether CGAL is used in the code or not
LIBCGAL=

## flags for building Tetgen, Deal.II and CGAL libraries
TETGEN_FLAGS=-DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY=../../lib 
DEALII_FLAGS=-DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DDEAL_II_WITH_NETCDF=OFF -DDEAL_II_STATIC_EXECUTABLE=ON -DDEAL_II_COMPONENT_DOCUMENTATION=OFF -DDEAL_II_COMPONENT_EXAMPLES=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../.. 
CGAL_FLAGS=-DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DBUILD_SHARED_LIBS=FALSE -DCMAKE_INSTALL_PREFIX=../.. 

## flag indicating in which machine the compilation happens
MACHINE=ubuntu

compiler=c++

main=Main.cpp

all: md_lib femocs_lib fecmd

fecmd: obj/fecmd.o
	${compiler} $< ${FEMOCS_OPT} ${FEMOCS_DOPT} ${FEMOCS_EXT_LIB} ${LIBPATH} ${LIB} -o $@

obj/fecmd.o: $(FEMOCS_PATH)/lib/libfemocs.a $(MD_PATH)/lib/libmd.a ${main}
	${compiler} -c ${main} ${FEMOCS_OPT} ${FEMOCS_DOPT} ${FEMOCS_EXT_HEAD} ${HEADPATH} -o $@

clean:
	rm -rf obj/*.o fecmd out/*

clean-all:
	rm -rf obj/*.o fecmd out/* && cd md && make clean && cd ../lib/femocs && make clean && cd GETELEC && make clean

md_lib:
	cd md && make md_lib

femocs_lib:
	cd lib/femocs/GETELEC && make && cd ../ && make lib
