
#ifdef _MAIN_

#else
#include <armadillo>
#include "../core/UAmoeba.hpp"
#ifdef _ALGEBRA_TOOLS_
#else
#include "../core/algebraTools.hpp"
#endif
#ifdef _AMOEBA_INIT_
#else
#include "../core/amoebaParam.cpp"
#endif

#endif

using namespace arma;


template<typename T>
void serialSolver(const AmoebaParam<T> amoebaParam)
{

	arma_rng::set_seed_random();

	vector<vec> newGuess;

	newGuess.resize(amoebaParam.nGridPoints);
	for(int i=0; i<amoebaParam.nGridPoints; ++i)
	{
		newGuess[i] = randu<vec>(amoebaParam.lieDimension);
	}

	UAmoeba* amoeba = new UAmoeba(amoebaParam.maxAmoebaIters,
																amoebaParam.nGridPoints,
																amoebaParam.precision,
																amoebaParam.matSize,
																amoebaParam.lieDimension,
																amoebaParam.basis);

	amoeba->solver(newGuess, amoebaParam.startBoundary, amoebaParam.endBoundary);
	amoeba->curvePrint();

}
