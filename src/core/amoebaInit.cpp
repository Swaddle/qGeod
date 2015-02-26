#define _AMOEBA_INIT_


#include "Pauli.hpp"

using namespace std;


template<typename S>
class AmoebaInit
{
public:
  int matSize;
  int numRows;

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
      matSize should be 2^n in SU(2^n)
      dimension = (2^n)^2 - 1 for SU(2^n)
      nGridPoints is number of starting guesses
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
    this->numRows = pauli->pauliBasisObject.size();

    cout << "Size of the lie algebra = " << numRows << "\n"
         << "Size of matrices = " << matSize << "\n";
  }
};
