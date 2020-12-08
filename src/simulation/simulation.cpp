#ifdef USE_MPI
#include <mpi.h>
#endif
#include <iostream>

#include "../state/array.hpp"
#include "../state/dynarray.hpp"
#include "../state/state.hpp"

//#include "../state/statefield.hpp"
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
  int num_proc=1, rank=0;
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Display MPI info
  int mpi_version,mpi_subversion;
  MPI_Get_version(&mpi_version, &mpi_subversion);
  if (rank == 0)
    std::cout << "Running with MPI version " << mpi_version << "."
              << mpi_subversion << " with " << num_proc << " processes"
              << std::endl;
#else
  std::cout << "Not using MPI\n";
#endif

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

  /*
  Array<double, 3> temp(0.0);
  Array<double, 3> temp2(1.0);

  std::function<Array<double, 3>(double, double, double)> f =
    [temp2] (double x, double y, double z) -> Array<double, 3> {
          return x*y*z*temp2;
        };

  temp += f(1, 1, 1);

  std::cout << temp[0] << std::endl;
  */

  State<UserSetup::nEq, UserSetup::maxOrder, UserSetup::nDim> S(mesh);
  MeshArray<Array<Array<double, UserSetup::nEq>, nDeg>> state(mesh->Nx, mesh->Ny, mesh->Nz);
  state = 0.0;

  // Set initial conditions
  InitialConditions<UserSetup::nEq> ic(cf);
  for (int i = mesh->nGhost; i < mesh->Nx - mesh->nGhost; i++)
    for (int j = mesh->nGhost; j < mesh->Ny - mesh->nGhost; j++)
      for (int k = mesh->nGhost; k < mesh->Nz - mesh->nGhost; k++)
        state[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i] = S.DoF(i, j, k, ic);

  /*
  RightHandSide<UserSetup::nEq, UserSetup::maxOrder, UserSetup::nDim> rhs(mesh);

  std::cout << "RightHandSide done" << std::endl;

  rhs.Calculate(state);
  //rhs.Calculate();
  */
  /*
  Array<double, UserSetup::nEq> B;
  B = 0.0;

  std::cout << "Hallo" << std::endl;

  auto f = [] (Array<double, UserSetup::nEq> A) { return A; };

  // Constructor, copy constructor, move constructor, move assignment  LEAK
  Array<double, UserSetup::nEq> result;
  result = f(B);

  std::cout << "Hallo " << result[0] << std::endl;
  */


  //Array<Array<double, 1>, 1> result;
  //result = 0.0;

  //std::cout << result[0][0] << std::endl;

  /*
  //TimeIntegrator<double, 3> ti;

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
  //delete[] state;
}

}
