/*
	* Implements leapfrog solving with multiple solvers
	*
	*
*/

#ifdef _MAIN_

#else
#include <armadillo>

#include "../core/UAmoeba.hpp"

#ifdef _ALGEBRA_TOOLS_
#else
#include "../core/algebraTools.hpp"
#endif

#ifdef _AMOEBA_INIT_
#else
#include "../core/amoebaParam.cpp"
#endif


#include "../core/translate.cpp"

using namespace algebraTools;
using namespace arma;

#endif

//cx_mat X0, cx_mat X1, int matSize, long maxAmoebaIters, long nGridPoints, double precision, long maxMainIters,

//size, rank MPI directives
//
template <typename T>
void segSolver(const AmoebaParam<T> amoebaParam, int rank, int size)
{

	cx_mat idMat = eye<T>(amoebaParam.matSize, amoebaParam.matSize);
	arma_rng::set_seed_random();


	// UAmoeba(long maxIters, int nGridPoints, double precision, int matSize, int lieDimension, vector<cx_mat> *inputBasis)

	cout << "Creating amoeba object \n";

	UAmoeba* amoeba = new UAmoeba(amoebaParam.maxAmoebaIters, amoebaParam.nGridPoints, amoebaParam.precision, amoebaParam.matSize, amoebaParam.lieDimension, amoebaParam.basis);

	cout << "Allocating MPI buffers \n";

	ArmadilloMPI* armaMPI = new ArmadilloMPI( amoebaParam.matSize, amoebaParam.matSize);
  vector<T> world;
  int bufferSize = 2 * size;

  if(rank == 0)
  {
   	//get X0
		cx_mat X0 = amoebaParam.startBoundary;
		cx_mat X1 = amoebaParam.endBoundary;

		cout << "Start Boundary = " << X0 << "\n"
				 << "Target Matrix = " << X1 << "\n";

    cx_mat kMid;
    invCayley<T>(X1, idMat, kMid);

   	world.push_back(X0);

		for(int i = 1; i < (bufferSize); ++i)
  	{
			kMid *= -0.5 * (i/static_cast<double>(bufferSize));
   		cayley<T>(kMid, idMat, X0);
  		world.push_back(X0);
  	}

    world.push_back(X1);
    cout << "Number of leap-frog points = "<< world.size() << "\n";
	}


	/*
		* int tag = 4 * rank;
		* need to send/recv 4 unique things
		* send and recieves have tag and tag+1 internally
		*
	*/

	int tag = 4 * rank;
	int offSet = 0;

	double localEnergy;
	double totalEnergy = 0;

	cx_mat bL(amoebaParam.matSize, amoebaParam.matSize);
	cx_mat bU(amoebaParam.matSize, amoebaParam.matSize);
	cx_mat newBound(amoebaParam.matSize, amoebaParam.matSize);


	vector<vec> newGuess;
	newGuess.resize(amoebaParam.nGridPoints);

	for(int i=0; i<amoebaParam.nGridPoints; ++i)
	{
		newGuess[i] = randu<vec>(amoebaParam.lieDimension);
	}
	/*
		* to keep track of where all the solvers are
		* each solver has to advance across the world
		* taking two boundary points
		* and outputing a new boundary
	*/

	int indLow;
	int indUp;
	int indMid;

	bufferSize = bufferSize + 1;

	int iters = 0;

	while(iters < amoebaParam.maxMainIters )
	{
		if(rank==0)
		{
			for(int src = 1; src < size; ++src)
			{
				indLow = (2 * src + offSet)%(bufferSize) ;
				indUp = (indLow + 2)%(bufferSize);
				if(indUp < indLow)
				{
					indLow = 0;
					indUp = indLow + 2;
				}

				//cout << "sending boundaries to " << src << endl;
				//cout << "sending " << indLow << " and " << indUp << endl;

				if(size > 0)
				{
					armaMPI->matDestroySend(world[indLow], src, 4 * src );
					armaMPI->matDestroySend(world[indUp], src, (4 * src) + 2);
				}
			}

			indLow = offSet%bufferSize;
			indUp = (indLow + 2)%(bufferSize);
			if(indUp < indLow)
			{
				indLow = 0;
				indUp = indLow + 2;
			}
			indMid = indLow + 1;
			//cout << "solving on " << indLow << " and " << indUp << endl;

			amoeba->curveSeeder(newGuess, amoebaParam.nGridPoints, world[indUp]);
			amoeba->solver(newGuess, world[indLow], world[indUp]);
			amoeba->newBoundary(world[indMid]);
			localEnergy = amoeba->getEnergy();

			for(int src = 1; src < size; ++src)
			{
				indLow = (2 * src + offSet)%(bufferSize) ;
				indUp = (indLow + 2)%(bufferSize);
				if(indUp < indLow)
				{
					indLow = 0;
				}
				indMid = indLow + 1;

				world[indMid] = armaMPI->matConstructRecv(src, 4 * src);
				//cout << "recieved updated boundary from" << src << endl;
			}

			offSet +=1;
			offSet = (offSet)%(bufferSize);
			//cout << "offSet is now " << offSet << endl;
		}
		else
		{
			bL = armaMPI->matConstructRecv(0, tag);
			bU = armaMPI->matConstructRecv(0, tag + 2);
			amoeba->curveSeeder(newGuess, amoebaParam.nGridPoints, bU);

			//cout << "recieved boundary on thread " << rank << endl;

			amoeba->solver(newGuess, bL, bU);
			amoeba->newBoundary(newBound);
			armaMPI->matDestroySend(newBound, 0, tag);
			localEnergy = amoeba->getEnergy();
		}

		MPI_Reduce(&localEnergy, &totalEnergy, size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if(rank == 0)
		{
			ofstream energyFile;
			energyFile.open("eOut.txt", ios::app);
			energyFile << totalEnergy << endl;
			cout << "iters" << iters << endl;
			cout << "=======================" << endl;
			totalEnergy = 0 ;

		}


		++iters;
	}

	if(rank == 0)
	{
		for(int i = 0; i < (2 * size - 1); ++i)
		{
			cout << "boundary is now" << world[i] << world[i+2];
			amoeba->curveSeeder(newGuess, amoebaParam.nGridPoints, world[i + 2]);
			amoeba->solver(newGuess, world[i], world[i + 2]);
			amoeba->curvePrint();
			amoeba->newBoundary(world[i+1]);

		}
	}

	delete armaMPI;
}
