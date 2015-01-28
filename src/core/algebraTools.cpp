
template <typename T>
void algebraTools::cayley(T &A, T &B, T &C)
{
  C = inv(A - B) * (A + B);
}

template <typename T>
void algebraTools::invCayley(T &A, T& B, T&C)
{
  C = (A - B) * (inv(A + B));
}

template <typename T>
void algebraTools::innerProd(T &A, T &B , double &scalar)
{
  scalar = real(trace(A.t() * B));
}

template <typename T>
void algebraTools::traceProd(T &A, T &B, double &scalar)
{
  scalar = real(trace(A * B));
}

template <typename T>
void lieBracket(T &A, T &B, T &C)
{
  C = A * B - B * A;
}
