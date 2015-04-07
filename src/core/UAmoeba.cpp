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

        for(int r = 0; r < lieDimension; ++r)
        {
            //kNew += kVector(r) * pauliBasis.pauliBasisObject[r];
            kNew += kVector(r) * (*basis)[r];

        }

        kNew *= -0.5 * h;
        cayley<cx_mat>(kNew, idMat, xNew);
        xNew *= xOld;
        xOld = xNew;
        kNew *= 0.0;
    }

    return xOld;
}

inline double UAmoeba::objectFunc(cx_mat A)
{
    return norm(A - endBoundary, 2);
}

inline double UAmoeba::cost(int index)
{
    if(index < 3)
    {
        return 10.0;
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


void UAmoeba::lieFunction(vec &vector)
{
    cx_mat k1 = zeros<cx_mat>(matSize, matSize);
    cx_mat k2 = zeros<cx_mat>(matSize, matSize);
    cx_mat k3;

    for(int t = 0; t < lieDimension; ++t)
    {
      k1 += vector(t) * (*basis)[t];
      k2 += vector(t) * cost(t) * (*basis)[t];

      //  k1 += vector(t) * pauliBasis.pauliBasisObject[t];
      //  k2 += vector(t) * cost(t) * pauliBasis.pauliBasisObject[t];
    }

    // computing [I k(t), k(t)]
    lieBracket<cx_mat>(k2, k1, k3);
    // return to vector form

    for(int s = 0; s < lieDimension-1; ++s)
    {
        // pauli basis is  0.5 i sigma_1, 0.5 i sigma_2 ....

        traceProd<cx_mat>(k3 , (*basis)[s] , vector(s));

      //  traceProd<cx_mat>(k3 , pauliBasis.pauliBasisObject[s], vector(s));
        vector(s) *= invCost(s);
    }

    vector *= -2 * h ;

    //note we need to scale this if not in SU(4)

    //    traceProd<cx_mat>(k3 , pauliBasis.pauliBasisObject[lieDimension-1], vector(lieDimension-1));


    traceProd<cx_mat>(k3 , (*basis)[lieDimension-1], vector(lieDimension-1));
    vector(lieDimension-1) *= h * (-1.0) * invCost(lieDimension-1);
}

void UAmoeba::amoebaRestart()
{
    //updates all but the best
    vec basisVec;

    for(int r = 1; r < lieDimension; ++r)
    {
      basisVec = zeros<vec>(lieDimension);
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
    vec guessVector = zeros<vec>(lieDimension);
    cx_mat K;
    invCayley<cx_mat>(eB, idMat, K);
    K *= 0.5;

    for(int i = 0; i < (lieDimension-1); ++i)
    {
      traceProd<cx_mat>(K , (*basis)[i], guessVector(i));
    }

    guessVector *= -2.0;

    traceProd<cx_mat>(K , (*basis)[lieDimension-1], guessVector(lieDimension-1) );

    for(int i = 0; i < nGridPoints; ++i)
    {
        newGuess[i] = guessVector + randu<vec>(lieDimension);
    }
}

/*
    * find the energy between boundary and end of guess curve
    * by joining with polynomial curve
*/

double UAmoeba::energyExtra(cx_mat A, cx_mat B)
{

    cx_mat W0;
    cx_mat tempMat = A.t() * B;

    invCayley<cx_mat>(tempMat , idMat, W0);
    cx_mat Z0 = A;
    cx_mat Z1(matSize, matSize);

    double localEnergy = 0.0;
    double temp;

    for(int r = 1; r < static_cast<double>(0.1 * gridSize); ++r)
    {
        W0 *= -0.5 * r * h * W0;
        cayley<cx_mat>(W0, idMat, Z0);
        Z0 *= Z0;

        Z1 = Z0 * Z0.t() - idMat;
        for(int p = 0; p < lieDimension; ++p)
        {
            innerProd<cx_mat>(Z1 , (*basis)[p], temp);
            localEnergy += pow(temp, 2);
        }
    }
    return localEnergy;

}

/*
  Evaluate the energy integral

  \int_{0}^{1} K(t) g K(t) dt

*/

void UAmoeba::amoebaEnergy()
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

        for(int r = 0; r < lieDimension; ++r)
        {
          kOld += kVectorOld(r) * (*basis)[r];
          kNew += kVector(r) * (*basis)[r];
          innerProd<cx_mat>(kNew, kOld, temp);
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
        for(int i = 0; i < (lieDimension - 1); ++i )
        {
            kFile << kVector(i) << ",";
        }
        kFile << kVector(lieDimension-1) << endl;

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

        for(int r = 0; r < lieDimension; ++r)
        {
            kNew += kVector(r) * (*basis)[r];
        }

        //xFile << real(trace(kNew.t() * kNew)) << endl;
        kNew *= -0.5 * h;
        cayley<cx_mat>(kNew, idMat, xNew);
        xNew *= xOld;
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
  for(int i = 0; i < (lieDimension-1); ++i)
  {
    traceProd<cx_mat>(B , (*basis)[i], v(i));
  }

  v *= -2.0;

  traceProd<cx_mat>(B , (*basis)[lieDimension-1], v(lieDimension-1));
  v(lieDimension-1) *= -1.0;

  for(int i = 0; i < (lieDimension - 1); ++i )
  {
      kFile << v(i) << ",";
  }
  kFile << v(lieDimension-1) << endl;

  cout << A * B ;
}
