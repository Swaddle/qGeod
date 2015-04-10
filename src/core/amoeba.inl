
/*
	* Implementation of amoeba
	*
*/

//=============================================================
/*
	* constructor
	* and destructor code
*/
//=============================================================

template<typename inType, typename outType>
Amoeba<inType, outType>::Amoeba(long maxIters, int dimension, double precision)
{
	//dimension is size of simplex (dimension of space + 1)
	this->dimension = dimension;
	this->maxIters = maxIters;
	this->precision = precision;

	//coefficients scale with dimension

	//reflect
	this->alpha = 1.0;

	//expansion
	this->beta = 1.0 +(2/static_cast<double>(dimension));

	//contraction
	this->gamma = 0.75 - (1/ 2 * static_cast<double>(dimension));

	//reduction
	this->delta = 1.0 - (1/static_cast<double>(dimension));

	wSimplex.resize(dimension);
	curveFuncResults.resize(dimension);
	this->last = dimension - 1;
}

template <typename inType, typename outType>
Amoeba<inType, outType>::~Amoeba()
{
	//stl vectors
}

//=============================================================
/*
	* solver
	* inputs a vector of initial points
	* outputs correct starting conditions
*/
//=============================================================

template <typename inType, typename outType>
void Amoeba<inType, outType>::solver(vector<inType>& initialGuess, outType startBoundary, outType endBoundary)
{
	srand(time(0));
	long iters = 0;
	long restartIters = 0;
	//long nRestarts = 0;
	double scaleFac = (1 /static_cast<double>(last));
	int secLast = dimension - 2;
	double objRefRes;
	double objExpRes;
	double objInConRes;
	double objOutConRes;
	double objAverage = 1000000;
	double objAverageOld = 0;
	double best;

	this->startBoundary = startBoundary;
	this->endBoundary = endBoundary;

	for(int i = 0; i < (dimension); ++i)
	{
		wSimplex[i].vertex = &initialGuess[i];
    curveFuncResults[i] = curveFunc(*wSimplex[i].vertex);
		wSimplex[i].weight = objectFunc(curveFuncResults[i]);
	}

	//sort the vertex based on sorted object results
	//sort the object function results

	__gnu_parallel::sort(wSimplex.begin(), wSimplex.end(), byWeight());


	//midpoint is calculated from all but the worst point
	midpoint = *wSimplex[0].vertex;
	for(int i = 1; i < secLast; ++i )
	{
		midpoint += *wSimplex[i].vertex;
	}
  midpoint *= scaleFac;
	best = wSimplex[0].weight;
	//inType midpointOld = midpoint;

	while(best>precision && iters<maxIters)
	{

		#pragma omp sections
		{
			{
				amoebaReflect();
				objRefRes = objectFunc(curveFunc(vertexReflected));
			}
			#pragma omp section
			{
				amoebaExpand();
				objExpRes = objectFunc(curveFunc(vertexExpanded));
			}
			#pragma omp section
			{
				amoebaOutSideContract();
				objOutConRes = objectFunc(curveFunc(vertexOutSideContracted));
			}

			#pragma omp section
			{
				amoebaInSideContract();
				objInConRes = objectFunc(curveFunc(vertexInSideContracted));
			}
		}



		if ( (wSimplex[0].weight <= objRefRes) && (objRefRes < wSimplex[secLast].weight) )
		{
			//cout << "[Ref]";
			*wSimplex[last].vertex = vertexReflected;
		}
		else if( objRefRes < wSimplex[0].weight)
		{

			*wSimplex[last].vertex = vertexReflected;

			if(objExpRes < objRefRes)
			{
				*wSimplex[last].vertex = vertexExpanded;
				//cout << "[E]";
			}
		}
		else if ( wSimplex[secLast].weight <= objRefRes && objRefRes < wSimplex[last].weight)
		{
			if(objOutConRes <= objRefRes)
			{
				//cout << "[OutC]";
				*wSimplex[last].vertex = vertexOutSideContracted;
			}
			else
			{
				//cout << "[Re]";
				amoebaReduce();
			}
		}
		else
		{
			// do inside contraction

			if(objInConRes < wSimplex[last].weight)
			{
				//cout << "[InC]";
				*wSimplex[last].vertex = vertexInSideContracted;
			}
			else
			{
				amoebaReduce();
			}
		}

		curveFuncResults[last] = curveFunc(*wSimplex[last].vertex);
		wSimplex[last].weight = objectFunc(curveFuncResults[last]);
		__gnu_parallel::sort(wSimplex.begin(), wSimplex.end(), byWeight());

		//return to start
		//re-calc midpoint
		//midpointOld = midpoint;

		midpoint = *wSimplex[0].vertex;
		objAverage = wSimplex[0].weight;
		for(int i = 1; i < secLast; ++i )
		{
			midpoint += *wSimplex[i].vertex;
			objAverage += wSimplex[i].weight;
		}
		midpoint *= scaleFac;
		objAverage *= scaleFac;

		//cout << "[O]" << objAverage;

		if( abs(objAverage - objAverageOld) < 0.00000001)
		{
			++restartIters;

			if(restartIters > 15)
			{
				//cout << " == RESTART == " << endl;
				amoebaRestart();

				midpoint = *wSimplex[0].vertex;
				objAverage = wSimplex[0].weight;
				for(int s = 1; s < last; ++s )
				{
					wSimplex[s].weight = objectFunc(curveFunc(*wSimplex[s].vertex));
					midpoint += *wSimplex[s].vertex;
					objAverage += wSimplex[s].weight;
				}
				midpoint *= scaleFac;
				objAverage *= scaleFac;
				std::sort(wSimplex.begin(), wSimplex.end(), byWeight());
				restartIters = 0;
			}
		}

		best = wSimplex[0].weight;
		objAverageOld = objAverage;
		++iters;
	}

	/*
	find energy of best vertex
	*/
	amoebaEnergy();


	cout << "====================" << endl;
	cout << " found in " << iters << " steps " << endl;
	cout << "optimal starting guess = " << endl;
	cout << *wSimplex[0].vertex << endl;
	cout << "vertex weight = " << endl;
	cout << wSimplex[0].weight << endl;
	cout << " X(1) = " << endl;
	cout << curveFunc(*wSimplex[0].vertex) << endl;
	cout << "====================" << endl;
}


//=============================================================
/*
	* functions that manipulate the simplex
	* dimension - 1 ---> c++ vectors are indexed at 0, we want
	* the worst vertex so the last in wSimplex list
*/
//=============================================================

template <typename inType, typename outType>
Amoeba<inType, outType>& Amoeba<inType, outType>::amoebaReflect()
{
	vertexReflected = ((1 + alpha) * midpoint) - (alpha * (*wSimplex[last].vertex));
	return *this;
}

template <typename inType, typename outType>
Amoeba<inType, outType>& Amoeba<inType, outType>::amoebaExpand()
{
	vertexExpanded = (( 1 - beta) * midpoint)  + (beta * vertexReflected);
	return *this;
}

template <typename inType, typename outType>
Amoeba<inType, outType>& Amoeba<inType, outType>::amoebaOutSideContract()
{
	vertexOutSideContracted = (( 1 - gamma) * midpoint) + (gamma * vertexReflected);
	return *this;
}

template <typename inType, typename outType>
Amoeba<inType, outType>& Amoeba<inType, outType>::amoebaInSideContract()
{
	vertexInSideContracted = (( 1 + gamma) * midpoint) - (gamma * vertexReflected);
	return *this;
}

template <typename inType, typename outType>
Amoeba<inType, outType>& Amoeba<inType, outType>::amoebaReduce()
{

	#pragma omp parallel for
		for(int i = 1; i < (dimension); ++i)
		{
			*wSimplex[i].vertex = ((1 - delta) * (*wSimplex[0].vertex))  + (delta * (*wSimplex[i].vertex));
		}
		
	return *this;
}

//=============================================================
/*
	* printer
	*
	*
*/
//=============================================================

template <typename inType, typename outType>
void Amoeba<inType, outType>::amoebaPrint(long index)
{

	//cout << "vertex value:" << wSimplex[index].vertex << endl;
	//cout << "curveFunc(vertex value) :" << curveFunc(wSimplex[index].vertex) << endl;
	//cout << "vertex weight :" << wSimplex[index].weight << endl;

}
