/*
  * qGeod
  *
*/

#define _MAIN_

#include <iostream>
#include <armadillo>


#include "../core/UAmoeba.hpp"


#include "../io/cxmatLoad.cpp"

#ifdef _ALGEBRA_TOOLS_
#else
#include "../core/algebraTools.hpp"
#endif

#ifdef _AMOEBA_INIT_
#else
#include "../core/amoebaParam.cpp"
#endif


#include "../core/translate.cpp"

#include "../solvers/segSolver.cpp"
#include "../solvers/serialSolver.cpp"
#include "../solvers/monte.cpp"

using namespace arma;

int main(int argc, char **argv)
{
	int rank;
	int size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	cx_mat startBoundary;
	cx_mat endBoundary;

	int maxIters;
	int precision;
	int matSize;
	int parallelFlag;

	if(argc==1)
	{
		cout<< "Usage : /qGeod <0/1/2> start.mat target.mat n maxAmoebaIters nGridPoints precision maxMainIters penalty\n"
				<< "=================================================================================================== \n"
				<< "n   				: exponent in SU(2^n)                                  \n"
				<< "precision   : working accuracy for qGeod                           \n"
				<< "maxIters    : maximum number of iterations in leap-frog            \n"
				<< "target.mat  : matrix containing desired unitary operation,         \n"
				<< "              in format described by                               \n"
			  << "              http://arma.sourceforge.net/docs.html#save_load_mat  \n"
				<< "nGridPoints : number of guess points in amoeba                     \n"
				<< "0/1/2       : run in serial 0, parallel 1, 2 monte carlo           \n"
				<< "=================================================================================================== \n";
	}
	else
	{

		parallelFlag = atoi(argv[1]);
		char *sBoundary = argv[2];
		char *eBoundary = argv[3];
		matSize = pow(2,atoi(argv[4]));

		cout << "Loading target matrix \n";
		cxmatLoad(endBoundary, matSize, eBoundary);
		cxmatLoad(startBoundary, matSize, sBoundary);

		cout << "Initialising solver parameters\n";
		AmoebaParam<cx_mat> amoebaParam;
		amoebaParam.getData(argc, argv);

		amoebaParam.startBoundary = startBoundary;
		amoebaParam.endBoundary = endBoundary;

		switch(parallelFlag)
		{
			case 0:
				cout << "Running serial solver\n";
				serialSolver<cx_mat>(amoebaParam);
				break;
			case 1:
				cout << "Running leap-frog solver\n";
				segSolver<cx_mat>(amoebaParam, rank, size);
				break;
			case 2:
				cout << "Running Monte Carlo solver\n";
				monteCarlo<cx_mat>(amoebaParam);
		}
	}


	MPI_Finalize();
	return 0;
}
