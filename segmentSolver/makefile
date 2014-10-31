
#set flags and compiler



CXX      	= mpic++ 
CXXFLAGS	= -std=c++11 -g -O2 -framework Accelerate -larmadillo 
#-O3 -Wall -I$(ARMAHOME)/include -L$(ARMAHOME)/lib -fopenmp -m64 -I$(INTEL)/include -L$(INTEL)/lib/intel64 $(INTELFLAGS) 
SRC 	= 
SOURCES  	=  ../solvers/segSolver.cpp
EXE = segSolver
#ARMAHOME =/ivec/sles11sp2/devel/gcc/4.4.7/armadillo/4.400.1
#LAPACKHOME = /ivec/sles11sp2/devel/gcc/4.4.7/lapack/3.4.2/lib
#INTEL = /ivec/intel/12.1.7/mkl
#INTELFLAGS = -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lmkl_blacs_intelmpi_lp64 -ldl -lpthread -lm
#MPIHOME = /opt/sgi/mpt/mpt-2.07/

all:
	$(CXX) $(SOURCES) $(CXXFLAGS) -o $(EXE) 









