#define _AMOEBA_INIT_


#include "Pauli.hpp"

using std::stod;
using std::atol;

template<typename S>
class AmoebaInit
{
public:
  int matSize;

  long maxAmoebaIters;
  long nGridPoints;
  long maxMainIters;
  std::vector<cx_mat> *basis;

  double precision;

  S startBoundary;
  S endBoundary;

  AmoebaInit()
  {
    //
  }

  void getData(int argc, char **argv)
  {

    /*
      matSize should be 2^n
      dimension = (2^n)^2 - 1 for SU(2^n)
    */
    this->matSize = atoi(argv[2]);


    this->maxAmoebaIters = atol(argv[3]);
    this->nGridPoints = atol(argv[4]);
    this->precision = stod(argv[5]);
    this->maxMainIters = stod(argv[6]);

    /*
      Using Pauli matrices which are 2^n x 2^n
      complex
      skew hermitian matrices
    */

    Pauli* pauli = new Pauli(matSize);

    this->basis = &pauli->pauliBasisObject;

  }
};
