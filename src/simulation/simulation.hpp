#ifndef DG_SIMULATION_HPP
#define DG_SIMULATION_HPP

#include "../state/array.hpp"
#include "../flux/flux.hpp"

namespace DGHydro {
  class Simulation {
  public:
    Simulation(char *fileName);
    ~Simulation();

  private:
    // Number of degrees of freedom
    const static int nDeg =
      (UserSetup::nDim == 1)*UserSetup::maxOrder +
      (UserSetup::nDim == 2)*(UserSetup::maxOrder + 1)*(UserSetup::maxOrder + 2)/2 +
      (UserSetup::nDim == 3)*(UserSetup::maxOrder + 1)*(UserSetup::maxOrder + 2)*(UserSetup::maxOrder + 3)/6;

    Array<Array<double, UserSetup::nEq>, nDeg> *state;

  };

}

#endif // DG_SIMULATION_HPP
