#ifndef DG_SIMULATION_HPP
#define DG_SIMULATION_HPP

#include "../array/array.hpp"
#include "../array/dynarray.hpp"
#include "../flux/flux.hpp"
#include "../state/state.hpp"

namespace DGHydro {
  using t_state = Array<double, UserSetup::nEq>;
  using t_state_deg = Array<Array<double, UserSetup::nEq>, UserSetup::nDeg>;

  class ConfigFile;
  class Mesh;

  class Simulation {
  public:
    Simulation(char *fileName);
    ~Simulation();

  private:
    ConfigFile *cf;
    Mesh *mesh;

    State<UserSetup::nEq, UserSetup::maxOrder, UserSetup::nDim> *state;

    DynArray<t_state_deg> *mesh_state;

    // Courant number
    double cfl;

    double CalcTimeStep();
  };

}

#endif // DG_SIMULATION_HPP
