#include <iostream>
#include <vector>
#include <armadillo>
#include <cmath>


using namespace std;
using namespace arma;


void inner(int depth,
          vector<int> & indices,
          vector<int> & maxIndex,
          vector<cx_mat>& Paulis,
          vector<cx_mat>& lieBasis
          )
{
  if (depth>0){
     for(int i = 0; i < maxIndex[depth-1]; ++i){
        indices[depth-1]=i;

        inner(depth-1, indices, maxIndex, Paulis, lieBasis);
     }
  }
  else
  {
    //do stuff with the index
      cout << "indices : ";

      cx_mat temp = Paulis[indices[0]];

      for(int r=0; r<indices.size(); ++r)
      {
        cout << indices[r] <<" ";
      }

      for(int r=1; r<indices.size(); ++r)
      {
        temp = kron(temp, Paulis[indices[r]]);
      }

      cout<<endl;

      lieBasis.push_back(temp);
   }
}

void nestedLoop(int n)
{

  vector<cx_mat> Paulis;
  Paulis.reserve(4);

  vector<cx_mat> lieBasis;

  Paulis[0] = eye<cx_mat>(2,2);
  Paulis[1] = { {0.0,0.0}, {1.0,0.0}, {1.0,0.0}, {0.0,0.0} };
  Paulis[1].reshape(2,2);
  Paulis[2] = { {0.0,0.0}, {0.0,1.0}, {0.0,-1.0}, {0.0,0.0} };
  Paulis[2].reshape(2,2);
  Paulis[3] = { {1.0 ,0.0}, {0.0,0.0}, {0.0,0.0}, {-1.0,0.0} };
  Paulis[3].reshape(2,2);

  vector<int>  indices(n,0);
  vector<int>  maxIndex;

  for(int i=0; i<n; ++i)
  {
    maxIndex.push_back(4);
  }

  inner(indices.size(),indices, maxIndex, Paulis, lieBasis);

  for(int i=0; i< lieBasis.size(); ++i)
  {

    cout << lieBasis[i];
    cout << endl;
  }

  cout << "size of lieBasis" << lieBasis.size();
}

int main(){
   nestedLoop(3);
}
