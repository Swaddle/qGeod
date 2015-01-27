/*
	* Implements leapfrog solving with multiple solvers
	*
	*
*/

#include "../core/UAmoeba.hpp"
#include "../core/translate.cpp"

int segSolver(cx_mat X0, cx_mat X1, int matSize, long maxAmoebaIters, long nGridPoints, double precision, long maxMainIters, int argc, char **argv)
{
	int rank;
	int size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	arma_rng::set_seed_random();


	UAmoeba* amoeba = new <matSize>(maxAmoebaIters, nGridPoints, precision);
  ArmadilloMPI* armaMPI = new ArmadilloMPI( matSize, matSize);

  vector<cx_mat> world;
  int bufferSize = 2 * size;

  if(rank == 0)
  {

   	//get X0


    cx_mat kMid;
    kMid = amoeba->invCayley(X1);

   	world.push_back(X0);
    for(int i = 1; i < (bufferSize); ++i)
    {
    	X0 = X0 * amoeba->cayley( -0.5 * (i/static_cast<double>(bufferSize)) * kMid);
    	world.push_back(X0);
    }
    world.push_back(X1);
    cout << world.size();
	}


	int iters = 0;
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

	cx_mat bL(4,4);
	cx_mat bU(4,4);
	cx_mat newBound(4,4);
	vector<vec> newGuess;
	newGuess.resize(nGridPoints);

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

	while(iters < maxMainIters )
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

				armaMPI->matDestroySend(world[indLow], src, 4 * src );
				armaMPI->matDestroySend(world[indUp], src, (4 * src) + 2);
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

			amoeba->curveSeeder(newGuess, nGridPoints, world[indUp]);
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
			amoeba->curveSeeder(newGuess, nGridPoints, bU);

			//cout << "recieved boundary on thread " << rank << endl;

			amoeba->solver(newGuess , bL, bU);
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
			amoeba->curveSeeder(newGuess, nGridPoints, world[i + 2]);
			amoeba->solver(newGuess, world[i], world[i + 2]);
			amoeba->curvePrint();
			amoeba->newBoundary(world[i+1]);

		}
	}

	delete armaMPI;
  MPI_Finalize();
  return 0;
}
