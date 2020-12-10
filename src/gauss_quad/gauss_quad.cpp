//#include <cmath>
#include <iostream>
#include <boost/math/special_functions/legendre.hpp>

#include "gauss_quad.hpp"

namespace DGHydro {

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Constructor
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  GaussQuad::GaussQuad(int n) : n(n)
{
  if (n <= 0)
    throw std::runtime_error("Can not create Gaussian Quadrature with a number of point that is not positive");

  x = new double[n];
  w = new double[n];

  FindAbscissae();
  FindWeights();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Destructor
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GaussQuad::~GaussQuad()
{
  delete[] x;
  delete[] w;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void GaussQuad::FindAbscissae()
{
  double dx = 2.0/(double) n;
  double tol = 1.0e-12;

  for (int i = 0; i < n; i++) {
    // Divide up interval: each contains a single zero
    double a = -1.0 + i*dx;
    double b = -1.0 + (i + 1)*dx;

    x[i] = 0.0;

    // Find zero by bisection
    while (1) {
      double c = 0.5*(a + b);

      double f = boost::math::legendre_p(n, c);

      if (f == 0.0 || 0.5*(b - a) < tol) {
        x[i] = c;
        break;
      }

      if (sgn(f) == sgn(boost::math::legendre_p(n, a))) {
        a = c;
      } else {
        b = c;
      }
    }

  }

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void GaussQuad::FindWeights()
{
  for (int i = 0; i < n; i++) {
    // Derivative of Legendre polynomial
    double p_prime = n*(x[i]*boost::math::legendre_p(n, x[i]) -
                        boost::math::legendre_p(n - 1, x[i]))/(x[i]*x[i]-1.0);
    // Gaussian weight
    w[i] = 2.0/(1 - x[i]*x[i])/p_prime/p_prime;
  }
}

}
