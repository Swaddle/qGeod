
#define _WEIGHTED_MAT_

#include <armadillo>

using arma::cx_mat;

struct weightedMat
{
  cx_mat lieMat;
  double cost;
};
