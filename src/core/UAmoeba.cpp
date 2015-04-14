/*
    * member functions for U()
    *
    *
*/

cx_mat UAmoeba::curveFunc(vec kVector)
{
    vec test = zeros<vec>(lieDimension);
    test(0) = 1;
    cout << "basis " << (*basis)[0].lieMat;
    cout << "Test " << matrixExp(test, 1.0 );


    //cx_mat kNew = zeros<cx_mat>(matSize, matSize);
    //cx_mat xNew = zeros<cx_mat>(matSize, matSize);
    cx_mat Xt = startBoundary;
    vec kVectorOld = kVector;

    for(int t = 0; t < gridSize; ++t)
    {

        Xt = Xt * matrixExp(kVectorOld,h);

        lieFunction(kVectorOld);
        kVector += h * kVectorOld;

        kVectorOld = kVector;

        /*for(int r = 0; r < lieDimension; ++r)
        //{
            //kNew += kVector(r) * pauliBasis.pauliBasisObject[r];
            //kNew += kVector(r) * (*basis)[r];

        //}

        //kNew *= -0.5 * h;
        //cayley<cx_mat>(kNew, idMat, xNew);
        */
        //xNew *= kNew *Xt;

        //Xt = xNew;
        //kNew *= 0.0;
    }
    return Xt;
}

inline double UAmoeba::objectFunc(cx_mat X1)
{
  return norm(X1 - endBoundary, 2);
}

inline double UAmoeba::cost(int index)
{
  return (*basis)[index].cost;
}

inline double UAmoeba::invCost(int index)
{
  return 1.0/((*basis)[index].cost);
}


void UAmoeba::lieFunction(vec &kVector)
{
    cx_mat k1 = kVector(0) * ((*basis)[0].lieMat);
    cx_mat k2 = kVector(0) * cost(0) * ((*basis)[0].lieMat);
    cx_mat k3;


    for(int r = 1; r < lieDimension; ++r)
    {
      k1 += kVector(r) * ((*basis)[r].lieMat);
      k2 += kVector(r) * cost(r) * ((*basis)[r].lieMat);
    }

    // computing k3 = [g k(t), k(t)]
    lieBracket<cx_mat>(k2, k1, k3);


    // return to vector form
    for(int s = 0; s < lieDimension; ++s)
    {
        traceProd<cx_mat>(k3 , (*basis)[s].lieMat , kVector(s));
        kVector(s) *= invCost(s);
    }

    //the trace gives you -8k_i , want -k_i cause right invariant metric has g k'(t) = -[g k(t),k(t)]
    kVector *= (1.0/matSize);
}

cx_mat UAmoeba::matrixExp(vec &kVector, double scalar)
{

  double kNorm = abs(norm(scalar * kVector,2));
  cx_mat cosPart = idMat * cos(kNorm);

  cx_mat sinPart = scalar * kVector(0) * (*basis)[0].lieMat;
  for(int r = 1; r<lieDimension; ++r)
  {
    sinPart += scalar * kVector(r) * (*basis)[r].lieMat;
  }
  sinPart *= (sin(kNorm)/kNorm);

  return (cosPart + sinPart);
}


void UAmoeba::amoebaRestart()
{
    //updates all but the best
    vec randVec;

    for(int r = 1; r < lieDimension; ++r)
    {
      randVec = zeros<vec>(lieDimension);
      randVec(r) = ((double) rand() / RAND_MAX);
      *wSimplex[r].vertex = (*wSimplex[0].vertex) + randVec;
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

    for(int i = 0; i < (lieDimension); ++i)
    {
      traceProd<cx_mat>(K , ((*basis)[i]).lieMat, guessVector(i));
    }

    guessVector *= (1.0/(matSize));

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
            innerProd<cx_mat>(Z1 , (*basis)[p].lieMat, temp);
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
          kOld += kVectorOld(r) * (*basis)[r].lieMat;
          kNew += kVector(r) * (*basis)[r].lieMat;
          innerProd<cx_mat>(kNew, kOld, temp);
          localE += pow(temp, 2);
          //cout << localE << endl;
        }

        kVectorOld = kVector;

        kNew *= 0.0;
        kOld *= 0.0;
    }

    //localE += energyExtra( Xt, endBoundary);

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
    kFile.open("kOut.csv",ios::app);
    xFile.open("xOut.csv",ios::app);

    cx_mat kNew = zeros<cx_mat>(matSize, matSize);
    cx_mat Xt = startBoundary;

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
              xFile << real(Xt(r,s)) << "," << imag(Xt(r,s)) << ",";
            }//
        }
        xFile << endl;

        Xt = Xt * matrixExp(kVectorOld, h);
        lieFunction(kVectorOld);
        kVector += kVectorOld;
        kVectorOld = kVector;
    }

    kFile.close();
    xFile.close();
}

void UAmoeba::curveCorrect(cx_mat &A, vec &v)
{
  cx_mat B =  0.5 * A.t() * endBoundary;
  for(int i = 0; i < (lieDimension-1); ++i)
  {
    traceProd<cx_mat>(B , (*basis)[i].lieMat, v(i));
  }

  v *= -2.0;

  traceProd<cx_mat>(B , (*basis)[lieDimension-1].lieMat, v(lieDimension-1));
  v(lieDimension-1) *= -1.0;

  for(int i = 0; i < (lieDimension - 1); ++i )
  {
      kFile << v(i) << ",";
  }
  kFile << v(lieDimension-1) << endl;

  cout << A * B ;
}
