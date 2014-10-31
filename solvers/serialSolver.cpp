#include "../src/UAmoeba.hpp"

int main(int argc, char **argv)
{

	long maxIters = atol(argv[1]);
	long nGridPoints = atoi(argv[2]);
	double precision = stod(argv[3]);
	arma_rng::set_seed_random();

	UAmoeba* amoeba = new UAmoeba(maxIters, nGridPoints, precision);

	cx_mat X0;
	X0 << cx_double(1.0, 0.0) << cx_double(0.0,0.0) << endr 
		 << cx_double(0.0,0.0) << cx_double(1.0, 0.0) << endr;

	cx_mat X1;
	X1 << cx_double(0.0, 0.0) << cx_double(0.0, 1.0 ) << endr
	    << cx_double(0.0, 1.0 ) << cx_double(0.0, 0.0) << endr;


	vector<vec> newGuess;
	for(int s = 0 ; s < nGridPoints; ++s)
	{
		newGuess.push_back(randu<vec>(4));
	}

	amoeba->solver(newGuess, X0, X1);
	amoeba->curvePrint();

}
