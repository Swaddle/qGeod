
/*
	* Generates the basis { i 0.5 sigma_1 , i 0.5 sigma_2 ......, i 0.5 sigma_n }
	* Where sigma_i is a generalized pauli matrix
*/

#define _PAULI
#ifdef _PAULI

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>

using namespace arma;
using std::vector;

class Pauli
{
	public:

		Pauli(int n);
		~Pauli();

		vector<cx_mat> pauliBasisObject;
		//cx_mat lieBracket(cx_mat A, cx_mat B);

	private:

		//size of matrix - n x n matrices
		long n;
		//dimension of space
		long dimension;
		int kD(int r, int s);

		cx_double cmplxI = cx_double(0.0, 0.5);
		cx_mat E(int r, int s);

		//generating methods
		Pauli& sigmaTypeOne();
		Pauli& sigmaTypeTwo();
		Pauli& sigmaTypeThree();

};

#include "Pauli.cpp"
#endif
