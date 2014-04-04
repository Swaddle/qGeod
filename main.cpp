

#include <cmath>
#include <iostream>

using namespace std; 

/*

	 1) Enumerate k with a guesses for the boundary say k(0) = I,  k(1) = InverseMatExp(U)

	 2) Shoot from the left and enumerate x'(i) = k(i) x(i) with the guess for k(0), x(0) = I find xLeft

	 3) Now shoot from the right and enumerate x'(i) = k(i) x(i) with the values for k(1), x(1) = U 
	  	find xRight
	
	 4) take the matrix norms of xRight and xLeft 

	 5) Evaluate xLeft - xRight = M 

	 6) Take k(0) =  I - InverseMatExp(M)
	




*/

double controlFunc(double x)
{

	return (x + 1) * x ; 
}

double inverse(double x);
{

	return (1/x); 
}

int main()
{

	int meshSize = 10000; 

	double k[meshSize], kOld[meshSize],x[meshSize], adj[meshSize];

	//initial condition, lets guess a k[0] 

	k[0] = 1;
	x[0] = 1; 

	double xTarget = 10;

	for(int i=1; i<= meshSize; i++)
	{


		k[i] = kOld[i-1] - (1 / meshSize) * controlFunc(kOld[i-1]); 
		x[i] = x[i-1] - (1/ meshSize) * ( kOld[i-1] * x[i-1] ) ;

		kOld[i] = k[i];

	}	


	if( abs(xTarget - x[meshSize]) > 0.001 ) 
	{

		//integrate back along adj(C(x, k))
		
		for(int i = 1; i<= meshSize; i++)
		{

			adj[i] = inverse(k[i])* x[i]* k[i];
		
		}


		x[0] = adj[meshSize]; 
	}
	else 
	{

		cout << "Solution Found" << endl;
	}


}