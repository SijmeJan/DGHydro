#ifndef DG_INITIAL_HPP
#define DG_INITIAL_HPP

#include <cmath>
#include <vector>

namespace DGHydro {
  class ConfigFile;

  class InitialConditions {
  public:
    InitialConditions(ConfigFile *cf) {};
    ~InitialConditions() {};

    std::vector<double> operator()(double x, double y, double z) {
      return std::vector<double>(1, cos(2.0*M_PI*x));
    }
  private:
};

} // namespace DGHydro

#endif  // DG_INITIAL_HPP
