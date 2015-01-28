
#include <iostream>
#include <armadillo>

using namespace arma;
using std::cout;
using std::endl; 

void cxmatLoad(cx_mat& uMat, int &size, char* file)
{
	cx_mat uTest;
	int uCheck;
	uTest.load(file, arma_ascii);

	if(uMat.n_rows !=  uMat.n_cols)
	{
		cout << "[Error] : Input matrix is non square" << endl;
		exit(1);
	}
	else
	{
		size = uMat.n_rows;
	}

	uCheck = static_cast<int>(real(trace(uMat * uMat.t())));

	if(uCheck != size)
	{
		cout << "[Error] : Input matrix is non-unitary " << endl;
		exit(1);
	}
	else
	{
		uMat = uTest;
	}
}
