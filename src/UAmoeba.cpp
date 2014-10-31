/*
    * member functions for U()
    *
    *
*/

cx_mat UAmoeba::curveFunc(vec kVector)
{
    cx_mat kNew = zeros<cx_mat>(matSize, matSize);
    cx_mat xNew = zeros<cx_mat>(matSize, matSize);
    cx_mat xOld = startBoundary;
    vec kVectorOld = kVector;

    for(int t = 1 ; t <= gridSize; ++t)
    {
        lieFunction(kVectorOld);
        kVector += kVectorOld;
        kVectorOld = kVector;

        for(int r = 0; r < numRows; ++r)
        {
            kNew += kVector(r) * _pauliBasis.pauliBasisObject[r];
        }

        xNew = xOld * cayley(-0.5 * h * kNew);
        xOld = xNew;
        kNew *= 0.0;
    }

    return xOld;
}

inline double UAmoeba::objectFunc(cx_mat A)
{
    return norm(A - endBoundary, 2);
}

cx_mat UAmoeba::cayley(cx_mat A)
{
    return (inv(iDMat - A)) * (iDMat + A);
}

cx_mat UAmoeba::invCayley(cx_mat A)
{
    return (iDMat - A) * (inv(iDMat + A));
}

inline double UAmoeba::cost(int index)
{
    if(index < 3)
    {
        return 100.0;
    }
    else
    {
        return 1.0;
    }
}

inline double UAmoeba::invCost(int index)
{
    if(index < 3 )
    {
        return 0.001;
    }
    else
    {
        return 1.0;
    }
}

double UAmoeba::innerProd(cx_mat &A, cx_mat &B)
{
  return real(trace(A.t() * B));
}

double UAmoeba::traceProd(cx_mat &A, cx_mat &B)
{
  return real(trace(A * B));
}


void UAmoeba::lieFunction(vec &vector)
{
    cx_mat k1 = zeros<cx_mat>(matSize, matSize);
    cx_mat k2 = zeros<cx_mat>(matSize, matSize);
    cx_mat k3 = zeros<cx_mat>(matSize, matSize);

    for(int t = 0; t < numRows; ++t)
    {
        k1 += vector(t) * _pauliBasis.pauliBasisObject[t];
        k2 += vector(t) * cost(t) * _pauliBasis.pauliBasisObject[t];
    }

    // computing [I k(t), k(t)]
    k3 = (k2 * k1) - (k1 * k2);
    // return to vector form

    for(int s = 0; s < numRows-1; ++s)
    {
        // pauli basis is  0.5 i sigma_1, 0.5 i sigma_2 ....
        vector(s) = invCost(s) * traceProd(k3 , _pauliBasis.pauliBasisObject[s]);
    }

    vector *= -2 * h ;

    //note we need to scale this if not in SU(4)
    vector(numRows-1) = h * (-1.0) * invCost(numRows-1) * traceProd(k3 , _pauliBasis.pauliBasisObject[numRows-1]);
}

void UAmoeba::amoebaRestart()
{
    //updates all but the best
    vec basisVec;

    for(int r = 1; r < numRows; ++r)
    {
      basisVec = zeros<vec>(numRows);
      basisVec(r) = ((double) rand() / RAND_MAX);
      *wSimplex[r].vertex = 0.75 * (*wSimplex[0].vertex) + basisVec;
    }

}

void UAmoeba::newBoundary(cx_mat& newBound)
{
    gridSizeOld = gridSize;
    gridSize = halfGridSize;
    newBound = curveFunc(*wSimplex[0].vertex);
    gridSize = gridSizeOld;
}


void UAmoeba::curveSeeder(vector<vec> &newGuess, int nGridPoints, cx_mat eB)
{
    vec guessVector = zeros<vec>(numRows);
    cx_mat K = 0.5 * invCayley(eB);

    for(int i = 0; i < (numRows-1); ++i)
    {
        guessVector(i) = (-2.0) * traceProd(K , _pauliBasis.pauliBasisObject[i]);
    }
    guessVector(numRows-1) = (-1.0) * traceProd(K , _pauliBasis.pauliBasisObject[numRows-1]);

    for(int i = 0; i < nGridPoints; ++i)
    {
        newGuess[i] = guessVector + randu<vec>(numRows);
    }
}

/*
    * find the energy between boundary and end of guess curve
    * by joining with polynomial curve
*/
double UAmoeba::energyExtra(cx_mat A, cx_mat B)
{

    cx_mat W0 = invCayley(A.t() * B);
    cx_mat Z0 = A;
    cx_mat Z1(matSize, matSize);

    double localEnergy = 0.0;

    for(int r = 1; r < static_cast<double>(0.1 * gridSize); ++r)
    {
        Z0 = Z0 * cayley( -0.5 * W0 * r * h);
        Z1 = Z0 * Z0.t() - iDMat;
        for(int p = 0; p < numRows; ++p)
        {
            localEnergy += pow(innerProd(Z1 , _pauliBasis.pauliBasisObject[p] ), 2);
        }
    }
    return localEnergy;

}

void UAmoeba::amoebaEnergy()
{
    vec kVector = *wSimplex[0].vertex;
    vec kVectorOld = kVector;
    cx_mat kNew = zeros<cx_mat>(matSize, matSize);
    cx_mat kOld = zeros<cx_mat>(matSize, matSize);

    double localE = 0.0;

    for(int t = 1 ; t <= gridSize; ++t)
    {
        lieFunction(kVectorOld);
        kVector += h * kVectorOld;

        for(int r = 0; r < numRows; ++r)
        {
          kOld += kVectorOld(r) * _pauliBasis.pauliBasisObject[r];
          kNew += kVector(r) * _pauliBasis.pauliBasisObject[r];
          localE += pow(innerProd(kNew, kOld), 2);
          //cout << localE << endl;
        }

        kVectorOld = kVector;

        kNew *= 0.0;
        kOld *= 0.0;
    }

    //localE += energyExtra( xOld, endBoundary);

    this->globalEnergy = localE;

    if( globalEnergy < globalEnergyOld)
    {
        globalEnergyOld = globalEnergy;
    }
    //eFile.open("eOut.txt",ios::app);
    //eFile << globalEnergy << endl;

}

double UAmoeba::getEnergy()
{
    //cout << "amoeba energy" << endl;
    //cout << globalEnergyOld << endl;
    return globalEnergyOld;
}


void UAmoeba::curvePrint()
{
    kFile.open("kOut.txt",ios::app);
    xFile.open("xOut.txt",ios::app);

    cx_mat kNew = zeros<cx_mat>(matSize, matSize);
    cx_mat xNew = zeros<cx_mat>(matSize, matSize);

    cx_mat xOld = startBoundary;
    vec kVector = *wSimplex[0].vertex;
    vec kVectorOld = kVector;
    //


    for(int t = 1 ; t <= gridSize; ++t)
    {
        for(int i = 0; i < (numRows - 1); ++i )
        {
            kFile << kVector(i) << ",";
        }
        kFile << kVector(numRows-1) << endl;

        for(int s = 0; s < matSize; ++s)
        {
          for(int r = 0; r < matSize; ++r)
            {
              xFile << real(xOld(r,s)) << "," << imag(xOld(r,s)) << ",";
            }//
        }
        xFile << endl;

        lieFunction(kVectorOld);
        kVector += kVectorOld;
        kVectorOld = kVector;

        for(int r = 0; r < numRows; ++r)
        {
            kNew += kVector(r) * _pauliBasis.pauliBasisObject[r];
        }

        //xFile << real(trace(kNew.t() * kNew)) << endl;

        xNew = xOld * cayley(-0.5 * h * kNew);
        xOld = xNew;
        kNew *= 0.0;


    }

    curveCorrect(xOld, kVector);


    kFile.close();
    xFile.close();
}

void UAmoeba::curveCorrect(cx_mat &A, vec &v)
{
  cx_mat B =  0.5 * A.t() * endBoundary;
  for(int i = 0; i < (numRows-1); ++i)
  {
    v(i) = (-2.0) * traceProd(B , _pauliBasis.pauliBasisObject[i]);
  }
  v(numRows-1) = (-1.0) * traceProd(B , _pauliBasis.pauliBasisObject[numRows-1]);

  for(int i = 0; i < (numRows - 1); ++i )
  {
      kFile << v(i) << ",";
  }
  kFile << v(numRows-1) << endl;

  cout << A * B ;
}
