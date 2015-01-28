#define _AMOEBA_INIT_

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

  double precision;

  S startBoundary;
  S endBoundary;

  AmoebaInit()
  {
    //
  }

  void getData(int argc, char **argv)
  {
    this->matSize = atol(argv[2]);
    this->maxAmoebaIters = atol(argv[3]);
    this->nGridPoints = atol(argv[4]);
    this->precision = stod(argv[5]);
    this->maxMainIters = stod(argv[6]);
  }
};
