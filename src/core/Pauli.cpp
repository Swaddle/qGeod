


//generate n by n matrices
Pauli::Pauli(long n)
{
	this -> n = n;

	dimension = (n * n) - 1;
	pauliBasisObject.reserve(dimension);

	cout << " No. of Pauli = " <<  dimension << endl;

	sigmaTypeThree();
	sigmaTypeOne();
	sigmaTypeTwo();
}

//kronecker delta
long Pauli::kD(long r, long s)
{
	return ~(r^s);
}

cx_mat Pauli::E(int r, int s)
{
	cx_mat A = zeros<cx_mat>(n,n);
	A(r - 1, s - 1) = cx_double(1.0, 0.0);
	return A;
}

//generating functions
Pauli& Pauli::sigmaTypeOne()
{
	for(int r = 2; r <= n; ++r)
	{
		for(int s = 1 ; s < r; ++s)
		{
			pauliBasisObject.push_back(cmplxI * (E(r, s) + E(s, r)));
		}
	}
	return *this;
}

Pauli& Pauli::sigmaTypeTwo()
{
	for(int r = 2; r <= n; ++r)
	{
		for(int s = 1 ; s < r ; ++s)
		{
			pauliBasisObject.push_back(E(s, r) - E(r, s));
		}
	}
	return *this;
}

Pauli& Pauli::sigmaTypeThree()
{
	for(int r = 1; r < n ; ++r )
	{
		cx_mat A = zeros<cx_mat>(n,n);

		for(int s = 1 ; s <= r  ; ++s)
		{
			A += E(s, s);
		}

		A -= r * E(r + 1, r + 1);
		pauliBasisObject.push_back( cmplxI * sqrt(2.0 / static_cast<double>(r * (r + 1))) * A);
	}
	return *this;
}

Pauli::~Pauli()
{

}
