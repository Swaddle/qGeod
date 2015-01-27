/*
  * qGeod
  *
*/

#include <iostream>
#include "../core/UAmoeba.hpp"
#include "../io/cxmatLoad.cpp"
#include "../solvers/segSolver.cpp"

int main(int argc, char **argv)
{

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
    cx_mat startMat = targetMat.t() * targetMat;

    int matSize = atol(argv[2]);
    long maxAmoebaIters = atol(argv[3]);
    long nGridPoints = atol(argv[4]);
    double precision = stod(argv[5]);
    long maxMainIters = stod(argv[6]);

    segSolver(startMat ,targetMat,  matSize,  maxAmoebaIters,  nGridPoints,  precision,  maxMainIters, argc, argv);

	}
	// load file

	//serial solve

	//parallel solve

	//exit
}
