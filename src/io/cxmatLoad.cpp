
using namespace armadillo

void cxmatLoad(cx_mat& uMat, int &size)
{
	cx_mat uTest;
	int uCheck; 
	uTest.load("../input/UnitaryMatrix.mat", arma_ascii); 

	if(uMat.n_rows != size)
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
