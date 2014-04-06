
/*
Geodesic solver for curves in SU(2^n)
Copyright (C) 2014  Michael Swaddle

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include <cmath>
#include <iostream>

using namespace std;


class Geodesic{




	public:
		
		
		int meshSize;

		Matrix x = new Matrix();
		Matrix k = new Matrix();




		double controlHamiltonian(double x)
		{

			return (x + 1) * x ; 
		}		





		
	
	private: 



};




int main()
{

	int meshSize = 10000; 

	Geodesic geod = new Geodesic(meshSize);


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