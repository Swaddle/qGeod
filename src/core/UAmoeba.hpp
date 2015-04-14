

//#define _DEBUG
#define _U_AMOEBA_

#ifdef _U_AMOEBA_
#include "amoeba.hpp"

#ifdef _ALGEBRA_TOOLS_
#else
#include "algebraTools.hpp"
#endif

#ifdef _WEIGHTED_MAT_
#else
#include "weightedMat.hpp"
#endif

#ifdef _LIE_ALGEBRA_
#else
#include "lieAlgebra.cpp"
#endif

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <armadillo>

using namespace algebraTools;
using std::vector;
using std::ofstream;
using std::complex;

using namespace arma;

//dimension
class UAmoeba: public Amoeba<vec, cx_mat>
{
public:
    UAmoeba(long maxIters,
            int nGridPoints,
            double precision,
            int matSize,
            int lieDimension,
            vector<weightedMat> *lieBasis)
    : Amoeba<vec, cx_mat>( maxIters,  nGridPoints,  precision)
    {

        //int matSize = number of rows in matrix
        //int lieDimension = number of basis vectors in lie algebra
        //int nGridPoints = number of guess points for the amoeba routine

        /*
          to do replace with parameter object
        */
        this->matSize = matSize;
        this->basis = lieBasis;
        this->gridSize = 1000;
        this->gridSizeOld = 1000;
        this->halfGridSize = 500;
        this->h = (1 / static_cast<double>(gridSize));
        this->idMat = eye<cx_mat>(matSize,matSize);
        this->lieDimension = lieDimension;
        this->globalEnergyOld = 10000000;

    }

    ~UAmoeba();

    //Pauli pauliBasis;
    void curvePrint();
    void newBoundary(cx_mat& newBound);
    void curveSeeder(vector<vec> &newGuess, int nGridPoints, cx_mat sU);
    double getEnergy();

protected:
    vector<weightedMat> *basis;

    int matSize;
    int gridSize, halfGridSize, gridSizeOld;
    int lieDimension;

    cx_double imagI = cx_double(0.0, 1.0);

    double globalEnergyOld;
    double globalEnergy;
    double h;

    cx_mat idMat;

    ofstream kFile;
    //ofstream eFile;
    ofstream xFile;

    //function describing the geodesic equations
    virtual cx_mat curveFunc(vec kVector);
    virtual double objectFunc(cx_mat A);
    virtual void amoebaEnergy();
    virtual void amoebaRestart();

    cx_mat matrixExp(vec &K, double scalar);

    double cost(int index);
    double invCost(int index);

    double energyExtra(cx_mat A, cx_mat B);
    void lieFunction(vec &kVector);
    void curveCorrect(cx_mat &A, vec &v);
};

#include "UAmoeba.cpp"
#endif
