

//#define _DEBUG
#define _U_AMOEBA_

#ifdef _U_AMOEBA_
#include "amoeba.hpp"


#ifdef _ALGEBRA_TOOLS_
#else
#include "algebraTools.hpp"
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

    UAmoeba(long maxIters, int dimension, double precision, vector<cx_mat> *inputBasis)
    : Amoeba<vec, cx_mat>( maxIters,  dimension,  precision)
    {
        //Pauli pauliBasis = new Pauli(dimension);
        this->matSize = dimension;
        this-> basis = inputBasis;
        //number of steps in integrators
        this->gridSize = 1000;
        this->gridSizeOld = 1000;
        this->halfGridSize = 500;
        this->h = (1 / static_cast<double>(gridSize));
        this->idMat = eye<cx_mat>(matSize,matSize);
        this->numRows = dimension; 
        this->globalEnergyOld = 10000000;
    }


    //Pauli pauliBasis;
    void curvePrint();
    void newBoundary(cx_mat& newBound);
    void curveSeeder(vector<vec> &newGuess, int nGridPoints, cx_mat sU);
    double getEnergy();

protected:
    vector<cx_mat> *basis;

    int matSize;
    int gridSize, halfGridSize, gridSizeOld;
    int numRows;

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

    double cost(int index);
    double invCost(int index);

    double energyExtra(cx_mat A, cx_mat B);
    void lieFunction(vec &kVector);
    void curveCorrect(cx_mat &A, vec &v);
};

#include "UAmoeba.cpp"
#endif
