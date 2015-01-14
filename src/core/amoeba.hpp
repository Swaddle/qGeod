
/*
	* Nelder Mead Simplex
*/

#define _AMOEBA
#ifdef _AMOEBA
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

//inType space of amoeba
//outTpe optimisation space
template<typename inType, typename outType>
class Amoeba
{
public:

	//object function exists in optimisation spacce
	virtual double objectFunc(outType) = 0;
	//function to map simplex to outType
	virtual outType curveFunc(inType) = 0;
	virtual void amoebaEnergy() = 0;
	virtual void amoebaRestart() = 0;

	Amoeba<inType, outType>(long maxIters, int dimension, double precision);
	~Amoeba<inType, outType>();
	void solver(vector<inType>& guess, outType startBoundary, outType endBoundary);

protected:

	struct weightedSimplex
	{
		inType *vertex;
		double weight;
	};

	struct byWeight
	{
    	bool operator()(weightedSimplex const &a, weightedSimplex const &b)
    	{
      		return a.weight < b.weight;
    	}

	};

	//the simplex object
	std::vector<weightedSimplex> wSimplex;
	std::vector<outType> curveFuncResults;

	Amoeba& amoebaReflect();
	Amoeba& amoebaExpand();
	Amoeba& amoebaInSideContract();
	Amoeba& amoebaOutSideContract();
	Amoeba& amoebaReduce();

	void amoebaPrint(long index);
	void midPoint(inType& midpoint);

	//points in the simplex
	inType vertexReflected;
	inType vertexExpanded;
	inType vertexInSideContracted;
	inType vertexOutSideContracted;
	inType midpoint;

	//critical points
	outType endBoundary;
	outType startBoundary;

	//scaling factors
	int dimension, last;
	long maxIters;
	double precision;
	double alpha, beta, gamma, delta;

};

#include "amoeba.inl"
#endif
