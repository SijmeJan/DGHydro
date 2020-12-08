#ifndef DG_FLUX_HPP
#define DG_FLUX_HPP

#include "../state/array.hpp"

namespace DGHydro {

  class UserSetup {
  public:
    UserSetup() {};
    ~UserSetup() {};

    const static int nDim = 3;         // Number of space dimensions
    const static int maxOrder = 0;     // Maximum order of polynomials
    const static int nEq = 1;          // Number of equations
  };

  // Scalar advection
  template<int nEq>
  class Flux {
  public:
    Flux() {};
    ~Flux() {};

    Array<double, nEq> x(const Array<double, nEq>& U) { return U;}
    Array<double, nEq> y(const Array<double, nEq>& U) { return U*0.0;}
    Array<double, nEq> z(const Array<double, nEq>& U) { return U*0.0;}

    double max_wave_speed_x(const Array<double, nEq>& U) {return 1.0;}
    double max_wave_speed_y(const Array<double, nEq>& U) {return 0.0;}
    double max_wave_speed_z(const Array<double, nEq>& U) {return 0.0;}
  private:
};

} // namespace DGHydro

#endif  // DG_FLUX_HPP
