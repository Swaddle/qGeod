/*
    * member functions for U()
    *
    *
*/

template <int DIMENSION>
cx_mat UAmoeba<DIMENSION>::curveFunc(vec kVector)
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

        kNew *= -0.5 * h;
        cayley(kNew, idMat, xNew);
        xNew *= xOld;
        xOld = xNew;
        kNew *= 0.0;
    }

    return xOld;
}

template <int DIMENSION>
inline double UAmoeba<DIMENSION>::objectFunc(cx_mat A)
{
    return norm(A - endBoundary, 2);
}

template <int DIMENSION>
inline double UAmoeba<DIMENSION>::cost(int index)
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

template <int DIMENSION>
inline double UAmoeba<DIMENSION>::invCost(int index)
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


template <int DIMENSION>
void UAmoeba<DIMENSION>::lieFunction(vec &vector)
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
    lieBracket(k2, k1, k3);
    // return to vector form

    for(int s = 0; s < numRows-1; ++s)
    {
        // pauli basis is  0.5 i sigma_1, 0.5 i sigma_2 ....
        traceProd(k3 , _pauliBasis.pauliBasisObject[s], vector(s));
        vector(s) *= invCost(s);
    }

    vector *= -2 * h ;

    //note we need to scale this if not in SU(4)
    traceProd(k3 , _pauliBasis.pauliBasisObject[numRows-1], vector(numRows-1));
    vector(numRows-1) *= h * (-1.0) * invCost(numRows-1);
}

template <int DIMENSION>
void UAmoeba<DIMENSION>::amoebaRestart()
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

template <int DIMENSION>
void UAmoeba<DIMENSION>::newBoundary(cx_mat& newBound)
{
    gridSizeOld = gridSize;
    gridSize = halfGridSize;
    newBound = curveFunc(*wSimplex[0].vertex);
    gridSize = gridSizeOld;
}


template <int DIMENSION>
void UAmoeba<DIMENSION>::curveSeeder(vector<vec> &newGuess, int nGridPoints, cx_mat eB)
{
    vec guessVector = zeros<vec>(numRows);
    cx_mat K;
    invCayley(eB, idMat, K);
    K *= 0.5;

    for(int i = 0; i < (numRows-1); ++i)
    {
        traceProd(K , _pauliBasis.pauliBasisObject[i], guessVector(i));
    }

    guessVector *= -2.0;

    traceProd(K , _pauliBasis.pauliBasisObject[numRows-1], guessVector(numRows-1) );

    for(int i = 0; i < nGridPoints; ++i)
    {
        newGuess[i] = guessVector + randu<vec>(numRows);
    }
}

/*
    * find the energy between boundary and end of guess curve
    * by joining with polynomial curve
*/

template <int DIMENSION>
double UAmoeba<DIMENSION>::energyExtra(cx_mat A, cx_mat B)
{

    cx_mat W0;
    cx_mat tempMat = A.t() * B;

    invCayley(tempMat , idMat, W0);
    cx_mat Z0 = A;
    cx_mat Z1(matSize, matSize);

    double localEnergy = 0.0;
    double temp;

    for(int r = 1; r < static_cast<double>(0.1 * gridSize); ++r)
    {
        W0 *= -0.5 * r * h * W0;
        cayley(W0, idMat, Z0);
        Z0 *= Z0;

        Z1 = Z0 * Z0.t() - idMat;
        for(int p = 0; p < numRows; ++p)
        {
            innerProd(Z1 , _pauliBasis.pauliBasisObject[p], temp);
            localEnergy += pow(temp, 2);
        }
    }
    return localEnergy;

}

template <int DIMENSION>
void UAmoeba<DIMENSION>::amoebaEnergy()
{
    vec kVector = *wSimplex[0].vertex;
    vec kVectorOld = kVector;
    cx_mat kNew = zeros<cx_mat>(matSize, matSize);
    cx_mat kOld = zeros<cx_mat>(matSize, matSize);

    double localE = 0.0;
    double temp;

    for(int t = 1 ; t <= gridSize; ++t)
    {
        lieFunction(kVectorOld);
        kVector += h * kVectorOld;

        for(int r = 0; r < numRows; ++r)
        {
          kOld += kVectorOld(r) * _pauliBasis.pauliBasisObject[r];
          kNew += kVector(r) * _pauliBasis.pauliBasisObject[r];
          innerProd(kNew, kOld, temp);
          localE += pow(temp, 2);
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

template <int DIMENSION>
double UAmoeba<DIMENSION>::getEnergy()
{
    //cout << "amoeba energy" << endl;
    //cout << globalEnergyOld << endl;
    return globalEnergyOld;
}

template <int DIMENSION>
void UAmoeba<DIMENSION>::curvePrint()
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
        kNew *= -0.5 * h;
        cayley(kNew, idMat, xNew);
        xNew *= xOld;
        xOld = xNew;
        kNew *= 0.0;


    }

    curveCorrect(xOld, kVector);


    kFile.close();
    xFile.close();
}

template <int DIMENSION>
void UAmoeba<DIMENSION>::curveCorrect(cx_mat &A, vec &v)
{
  cx_mat B =  0.5 * A.t() * endBoundary;
  for(int i = 0; i < (numRows-1); ++i)
  {
    traceProd(B , _pauliBasis.pauliBasisObject[i], v(i));
  }

  v *= -2.0;

  traceProd(B , _pauliBasis.pauliBasisObject[numRows-1],v(numRows-1));
  v(numRows-1) *= -1.0;

  for(int i = 0; i < (numRows - 1); ++i )
  {
      kFile << v(i) << ",";
  }
  kFile << v(numRows-1) << endl;

  cout << A * B ;
}
