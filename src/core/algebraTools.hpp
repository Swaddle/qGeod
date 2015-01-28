#define _ALGEBRA_TOOLS_

#ifdef _ALGEBRA_TOOLS_

namespace algebraTools
{

  template <typename T>
  void cayley(T &A, T &B, T &C);

  template <typename T>
  void invCayley(T &A, T& B, T&C);

  template <typename T>
  void innerProd(T &A, T &B , double &scalar);

  template <typename T>
  void lieBracket(T &A, T &B, T &C);

  template <typename T>
  void traceProd(T &A, T &B, double &scalar);
}

#include "algebraTools.cpp"
#endif
