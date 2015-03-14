/*
  * qGeod
  *
*/

#include <iostream>
#include "../core/amoebaInit.cpp"
#include "../io/cxmatLoad.cpp"
#include "../solvers/segSolver.cpp"
#include "armadillo"


using namespace arma;

int main(int argc, char **argv)
{
	int rank;
	int size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	cx_mat targetMat;
	cx_mat initMat;

	int maxIters;
	int precision;
	int matSize;

	if(argc==1)
	{
		cout<< "Usage : /qGeod start.mat target.mat matSize maxAmoebaIters nGridPoints precision maxMainIters \n"
				<< "================================================================= \n"
				<< "precision  : working accuracy for qGeod     \n"
				<< "maxIters   : maximum number of iterations in leap-frog \n"
				<< "target.mat : matrix containing desired unitary operation, \n"
				<< "             in format described by \n"
			  << "						 http://arma.sourceforge.net/docs.html#save_load_mat \n"
				<< "================================================================= \n";
	}
	else
	{

	char *startMat = argv[1];
	char *stopMat = argv[2];

	matSize = atoi(argv[3]);

	cout << " Loading target matrix \n";


	cxmatLoad(targetMat, matSize, stopMat);
	cxmatLoad(initMat, matSize, startMat);

	cout << "Initialising solver parameters\n";

	AmoebaInit<cx_mat> amoebaParam;
	amoebaParam.getData(argc, argv);
	amoebaParam.startBoundary = initMat;
	amoebaParam.endBoundary = targetMat;

	cout << "Running leap-frog solver\n";


	segSolver<cx_mat>(amoebaParam, rank, size);

	}

	MPI_Finalize();
	return 0;
}
