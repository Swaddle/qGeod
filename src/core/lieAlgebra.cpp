
#define _LIE_ALGEBRA_


#include <iostream>
#include <vector>
#include <armadillo>
#include <cmath>


#ifdef _WEIGHTED_MAT_
#else
#include "weightedMat.hpp"
#endif

using namespace std;
using namespace arma;

class LieAlgebra
{

  public:

    LieAlgebra(int n, double penalty)
    {
      basisGenerator(n, penalty, lieBasis);
    }

    ~LieAlgebra();

    vector<weightedMat> lieBasis;

  private:

    cx_double imagI = cx_double(0.0, 1.0);

    void generatorLoop(int depth,
              vector<int> & indices,
              vector<int> & maxIndex,
              vector<cx_mat>& Paulis,
              vector<weightedMat>& lieBasis,
              double penalty
              )
    {
      if (depth>0){
         for(int i = 0; i < maxIndex[depth-1]; ++i){
            indices[depth-1] = i;
            generatorLoop(depth-1, indices, maxIndex, Paulis, lieBasis, penalty);
         }
      }
      else
      {
          // indices of tensor product
          //cout << "indices : ";
          int numNonZeros = 0;
          weightedMat temp;
          temp.lieMat = Paulis[indices[0]];

          for(int r = 0; r < indices.size(); ++r)
          {
            if(indices[r] != 0)
            {
              ++numNonZeros;
            }
            //cout << indices[r] <<" ";
          }

          if(numNonZeros < 3 )
          {
            temp.cost = 1;
          }
          else
          {
            temp.cost = penalty;
          }

          //make kronecker products
          for(int r = 1; r < indices.size(); ++r)
          {
            temp.lieMat = kron(temp.lieMat, Paulis[indices[r]]);
          }

          temp.lieMat *= imagI;

          lieBasis.push_back(temp);
       }
    }

    //n is 2^n in SU(2^n)
    void basisGenerator(int n, double penalty, vector<weightedMat>& lieBasis)
    {

      vector<cx_mat> Paulis;
      Paulis.reserve(4);

      Paulis[0] = eye<cx_mat>(2,2);
      Paulis[1] = { {0.0,0.0}, {1.0,0.0}, {1.0,0.0}, {0.0,0.0} };
      Paulis[2] = { {0.0,0.0}, {0.0,1.0}, {0.0,-1.0}, {0.0,0.0} };
      Paulis[3] = { {1.0,0.0}, {0.0,0.0}, {0.0,0.0}, {-1.0,0.0} };

      for(int s = 1; s < 4; ++s)
      {
        Paulis[s].reshape(2,2);
      }

      vector<int> indices(n,0);
      vector<int> maxIndex;

      for(int i=0; i<n; ++i)
      {
        maxIndex.push_back(4);
      }

      generatorLoop(indices.size(),indices, maxIndex, Paulis, lieBasis, penalty);

      //for SU remove the identity product I otimes I .... otimes I
      lieBasis.erase(lieBasis.begin());
    }
};
