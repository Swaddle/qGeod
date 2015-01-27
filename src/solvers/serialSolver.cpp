#include "../core/UAmoeba.hpp"

int serialSolver(int dimension)
{

	long maxIters = atol(argv[1]);
	long nGridPoints = atoi(argv[2]);
	double precision = stod(argv[3]);
	arma_rng::set_seed_random();

	UAmoeba<4>* amoeba = new UAmoeba<dimension>(maxIters, nGridPoints, precision);

	//get boundary points

	vector<vec> newGuess;
	for(int s = 0 ; s < nGridPoints; ++s)
	{
		newGuess.push_back(randu<vec>(4));
	}

	amoeba->solver(newGuess, X0, X1);
	amoeba->curvePrint();

}
