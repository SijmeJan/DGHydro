#ifndef DG_BASIS_HPP
#define DG_BASIS_HPP

#include <tuple>
#include <map>
#include <boost/math/special_functions/legendre.hpp>

namespace DGHydro {

  template<int nDim, int nDeg>
  class BasisFunctions {
  public:
    BasisFunctions() {
      total_number = 0;

      for (int i = 0; i <= nDeg; i++) {
        for (int j = 0; j <= nDeg*(nDim > 1); j++) {
          for (int k = 0; k <= nDeg*(nDim > 2); k++) {
            if (i + j + k <= nDeg) {
              basisMap.insert({total_number, std::make_tuple(i, j, k)});

              std::cout << "Added " << i << " " << j << " " << k << std::endl;
              total_number++;
            }
          }
        }
      }
    };

    ~BasisFunctions(void) {};

    double operator()(int i, double x, double y, double z) {
      int a = std::get<0>(basisMap[i]);
      int b = std::get<1>(basisMap[i]);
      int c = std::get<2>(basisMap[i]);

      return boost::math::legendre_p(a, x)*boost::math::legendre_p(b, y)*boost::math::legendre_p(c, z);
    }

    double x_derivative(int i, double x, double y, double z) {
      int a = std::get<0>(basisMap[i]);
      int b = std::get<1>(basisMap[i]);
      int c = std::get<2>(basisMap[i]);

      double p_prime = a*(x*boost::math::legendre_p(a, x) -
                          boost::math::legendre_p(a - 1, x))/(x*x - 1.0);

      return p_prime*boost::math::legendre_p(b, y)*boost::math::legendre_p(c, z);
    }

    double y_derivative(int i, double x, double y, double z) {
      int a = std::get<0>(basisMap[i]);
      int b = std::get<1>(basisMap[i]);
      int c = std::get<2>(basisMap[i]);

      double p_prime = b*(y*boost::math::legendre_p(b, y) -
                          boost::math::legendre_p(b - 1, y))/(y*y - 1.0);

      return boost::math::legendre_p(a, x)*p_prime*boost::math::legendre_p(c, z);
    }

    double z_derivative(int i, double x, double y, double z) {
      int a = std::get<0>(basisMap[i]);
      int b = std::get<1>(basisMap[i]);
      int c = std::get<2>(basisMap[i]);

      double p_prime = c*(z*boost::math::legendre_p(c, z) -
                          boost::math::legendre_p(c - 1, z))/(z*z - 1.0);

      return boost::math::legendre_p(a, x)*boost::math::legendre_p(b, y)*p_prime;
    }


    int total_number;
  private:
    std::map<int, std::tuple<int, int, int>> basisMap;


  };

} // namespace DGHydro

#endif  // DG_BASIS_HPP
