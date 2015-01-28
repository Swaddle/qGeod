/*
  * qGeod
  *
*/

#include <iostream>
#include "../core/amoebaInit.cpp"
#include "../io/cxmatLoad.cpp"
#include "../solvers/segSolver.cpp"


int main(int argc, char **argv)
{
	int rank;
	int size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	cx_mat targetMat;
	int maxIters;
	int precision;
	int matSize;

	if(argc==1)
	{
		cout<< "Usage : /qGeod <-s> <-p> target.mat  matSize, maxAmoebaIters,  nGridPoints,  precision, maxMainIters \n"
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

		char *filename = argv[1];
		cxmatLoad(targetMat, matSize, filename);

		AmoebaInit<cx_mat> amoebaParam;
		amoebaParam.getData(argc, argv);
		amoebaParam.startBoundary = eye<cx_mat>(matSize, matSize);
		amoebaParam.endBoundary = targetMat;

    segSolver<cx_mat>(amoebaParam, rank, size);

	}
	// load file

	//serial solve

	//parallel solve

	//exit
}
