#include <cstdlib>
#include <ctime>
#include <iostream>
#include <armadillo>
#include <random>
#include <functional>


using namespace arma;


class RandDouble
{
  public:
      RandDouble(double low, double high)
      :r(std::bind(std::uniform_real_distribution<>(low,high),std::default_random_engine())){}

      double operator()(){ return r(); }

  private:
      std::function<double()> r;
};


int randInt(int minInt, int  maxInt)
{
  return rand()%maxInt + minInt;
}


template <typename T>
void monteCarlo(const AmoebaParam<T> amoebaParam)
{

  srand(time(NULL));

  vec guess(63);
  for(int i=0; i<63; ++i)
  {
      if(i<37){
        guess(i) = 0.0001;
      }
      else
      {
        guess(i) = 0.000000000001;
      }
  }
  vec newGuess = guess;

  cx_mat X1;

  double tau = 6.28318530718;
  double tol = 0.000001;
  double objectRes = 100;
  double objectResOld;

  double change;

  int entry;
  int iters = 0;
  int maxIters = 10000;
  double scale = 1.0;

  RandDouble randD{-0.00001,0.00001};

  UAmoeba* amoeba = new UAmoeba(amoebaParam.maxAmoebaIters, amoebaParam.nGridPoints, amoebaParam.precision, amoebaParam.matSize, amoebaParam.lieDimension, amoebaParam.basis);
  amoeba->startBoundary = amoebaParam.startBoundary;
  amoeba->endBoundary = amoebaParam.endBoundary;

  cout << "start boundary : " << amoeba->startBoundary << endl;
  cout << "end boundary : " << amoeba->endBoundary << endl;

  X1 = amoeba->curveFunc(guess);
  objectResOld = amoeba->objectFunc(X1);

  cout << "initial weight : " << objectResOld << endl;

  while( iters < maxIters && objectRes > tol)
  {

    iters+=1;
    cout << iters << endl;
    //randomly select a entry
    entry = iters%38;
    //randomly vary the entry mod 2pi
    change = randD();
    newGuess(entry) = fmod((newGuess(entry) + change), tau);

    //recompute geodesic
    X1 = amoeba->curveFunc(newGuess);

    //compute norm
    objectRes = amoeba->objectFunc(X1);


    if(objectRes < objectResOld)
    {
      //new guess is better
      cout << "changed entry :" << entry << endl;
      cout <<"change : " << change << endl;
      cout << "new weight : " << objectRes << endl;
      guess = newGuess;
      objectResOld = objectRes;
    }
    //reset newGuess;
    newGuess = guess;
  }

  amoeba->curvePrint(guess);
  cout << "best guess : " << guess << endl;
}
