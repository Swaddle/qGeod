

//#define _DEBUG
#define _U_AMOEBA_

#ifdef _U_AMOEBA_

#include "amoeba.hpp"
#include "Pauli.hpp"
#include "algebraTools.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <complex>

using namespace std;
using namespace algebraTools;

using std::vector;

//dimension
template<int DIMENSION> class UAmoeba: public Amoeba<vec, cx_mat>
{
public:

    UAmoeba(long maxIters, int dimension, double precision)
    : Amoeba<vec, cx_mat>( maxIters,  dimension,  precision)
    , _pauliBasis(DIMENSION)
    {
        this->matSize = DIMENSION;
        //number of steps in integrators
        this->gridSize = 1000;
        this->gridSizeOld = 1000;
        this->halfGridSize = 500;
        this->h = (1 / static_cast<double>(gridSize));

        this->idMat = zeros<cx_mat>(matSize,matSize);
        for(int r = 0; r < matSize; ++r )
        {
          this->idMat(r,r) = cx_double(1.0,0.0);
        }
        _pauliBasis.pauliBasisObject.push_back(cx_double(0.0,0.5)*idMat);
        this->numRows = (int)_pauliBasis.pauliBasisObject.size();
        this->globalEnergyOld = 10000000;
    }

    Pauli _pauliBasis;
    void curvePrint();
    void newBoundary(cx_mat& newBound);
    void curveSeeder(vector<vec> &newGuess, int nGridPoints, cx_mat sU);
    double getEnergy();

protected:

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
