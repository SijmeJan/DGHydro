#ifndef DG_INITIAL_HPP
#define DG_INITIAL_HPP

#include <cmath>
//#include <vector>

#include "../state/array.hpp"

namespace DGHydro {
  class ConfigFile;

  template<int nEq>
  class InitialConditions {
  public:
    InitialConditions(ConfigFile *cf) {};
    ~InitialConditions() {};

    Array<double, nEq> operator()(double x, double y, double z) {
      return Array<double, nEq>(x);//(cos(2.0*M_PI*x));
    }
  private:
};

} // namespace DGHydro

#endif  // DG_INITIAL_HPP
