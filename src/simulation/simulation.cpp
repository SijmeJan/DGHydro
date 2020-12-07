#include <mpi.h>
#include <iostream>

#include "../state/array.hpp"
#include "../state/dynarray.hpp"

#include "../state/statefield.hpp"
#include "../state/basis.hpp"
#include "simulation.hpp"
#include "../config_file/config_file.hpp"
#include "../mesh/mesh.hpp"
#include "../gauss_quad/gauss_quad.hpp"
#include "../gauss_quad/integral.hpp"
#include "../initial/initial.hpp"
#include "../rhs/rhs.hpp"
#include "../flux/flux.hpp"
#include "time_integrator.hpp"
#include "../state/mesharray.hpp"

namespace DGHydro {

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Constructor
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Simulation::Simulation(char *fileName)
{
  /*
  // Defines number of space dimensions etc
  const UserSetup setup;

  // Number of degrees of freedom
  const int nDeg =
    (setup.nDim == 1)*setup.maxOrder +
    (setup.nDim == 2)*(setup.maxOrder + 1)*(setup.maxOrder + 2)/2 +
    (setup.nDim == 3)*(setup.maxOrder + 1)*(setup.maxOrder + 2)*(setup.maxOrder + 3)/6;
  */

  // Initialize fftw3
  //fftw_mpi_init();

  // Get total number of CPUs and rank of current CPU
  int num_proc, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Display MPI info
  int mpi_version,mpi_subversion;
  MPI_Get_version(&mpi_version, &mpi_subversion);
  if (rank == 0)
    std::cout << "Running with MPI version " << mpi_version << "."
              << mpi_subversion << " with " << num_proc << " processes"
              << std::endl;

  // Read configuration file
  ConfigFile *cf;
  try {
    cf = new ConfigFile(fileName, rank);
    cf->List();
  }
  catch (...) {
    std::cout << "Error reading configuration file!" << std::endl;
    throw;
  }

  // Build mesh
  Mesh *mesh;
  try {
    mesh = new Mesh(cf);
    mesh->Decompose(rank, num_proc);
  }
  catch (std::exception& e) {
    std::cout << e.what() << '\n';
    throw std::runtime_error("Could not create simulation");
  }

  //state = new StateField<1,1>(mesh->Nx, mesh->Ny, mesh->Nz);

  //InitialConditions ic = InitialConditions(cf);
  //for (int i = 0; i < mesh->Nx; i++) {
  //  state->get<0,0>(i) = ic(mesh->x[i], 0, 0)[0];
  //}

  //for (int i = 0; i < mesh->Nx; i++) {
  //  std::cout << state->get<0,0>(i) << " ";
  //}

  //GaussQuad *gq = new GaussQuad(5);
  //const int N = 1;
  //CubeIntegral<State<N>, 5> *ci = new CubeIntegral<State<N>, 5>;
  //State<N> s;

  //s = 2*(s + 1.0);

  //std::cout << ([s](double x) -> State<N> { return s + x; })(4.0)[0] << std::endl;

  //std::function<State<N>(double)> f = [&s] (double x) -> State<N> { return State<N>(s + 1.0);};

  //State<N> res = ci->vol1d(f);
  //std::cout << res[0] << std::endl;


  //BasisFunctions<1, 2> bf;

  //std::cout << bf.x_derivative(0, 0, 0, 0) << std::endl;

  //state = new Array<Array<double, UserSetup::nEq>, nDeg>[mesh->Nx*mesh->Ny*mesh->Nz];

  MeshArray<Array<Array<double, UserSetup::nEq>, nDeg>> state(mesh->Nx, mesh->Ny, mesh->Nz);

  RightHandSide<UserSetup::nEq, UserSetup::maxOrder, UserSetup::nDim> rhs(mesh);

  rhs.Calculate(state);

  std::cout << rhs.data[0][0][0][0] << std::endl;

  /*
  TimeIntegrator<double, 3> ti;

  double u = 1.0;
  for (int i = 0; i < 10; i++)
    ti.TakeStep(0.0, 0.1, u, [](double t, double u){ return u; });

  std::cout << u << std::endl;
  */

  delete cf;
  delete mesh;
}

Simulation::~Simulation()
{
  delete[] state;
}

}
