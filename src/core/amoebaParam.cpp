#define _AMOEBA_INIT_

#ifdef _LIE_ALGEBRA_
#else
#include "lieAlgebra.cpp"
#endif

using std::vector;
using std::cout;
using std::stod;
using std::atol;
using std::pow;

template<typename S>
class AmoebaParam
{
public:
  int matSize;
  int lieDimension;

  long maxAmoebaIters;
  long nGridPoints;
  long maxMainIters;
  std::vector<weightedMat> *basis;

  double precision;
  double penalty;

  S startBoundary;
  S endBoundary;

  AmoebaParam()
  {
    //
  }

  void getData(int argc, char **argv)
  {

    /*
      matSize should be 2^n in SU(2^n)
      dimension = (2^n)^2 - 1 for SU(2^n)
      nGridPoints is number of starting guesses for Nelder Mead
    */
    this->matSize = pow(2,atoi(argv[4]));
    this->maxAmoebaIters = atol(argv[5]);
    this->nGridPoints = atol(argv[6]);
    this->precision = stod(argv[7]);
    this->maxMainIters = stod(argv[8]);
    this->penalty = stod(argv[9]);

    /*
      Using Pauli matrices which are 2^n x 2^n
      complex
      skew hermitian matrices

      LieAlgebra(n, penalty)
    */

    LieAlgebra* lieAlgebra = new LieAlgebra(atoi(argv[4]), penalty);

    this->basis = &lieAlgebra->lieBasis;

    // SU(d) - dimension d^2
    this->lieDimension = lieAlgebra->lieBasis.size();

    cout << "Size of the lie algebra = " << lieDimension << "\n";
  }


  ~AmoebaParam()
  {
  }

};
