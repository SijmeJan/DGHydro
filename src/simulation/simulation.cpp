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
    std::cout << "Hallo\n";
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

  State<UserSetup::nEq, UserSetup::maxOrder, UserSetup::nDim>
    state(mesh);
  DynArray<Array<Array<double, UserSetup::nEq>, nDeg>>
    mesh_state(mesh->Nx*mesh->Ny*mesh->Nz);
  mesh_state = 0.0;

  // Set initial conditions
  InitialConditions<UserSetup::nEq> ic(cf);
  for (int i = mesh->nGhost; i < mesh->Nx - mesh->nGhost; i++)
    for (int j = mesh->nGhost; j < mesh->Ny - mesh->nGhost; j++)
      for (int k = mesh->nGhost; k < mesh->Nz - mesh->nGhost; k++)
        mesh_state[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i] =
          state.DoF(i, j, k, ic);

  RightHandSide<UserSetup::nEq, UserSetup::maxOrder, UserSetup::nDim> rhs(mesh);

  TimeIntegrator<DynArray<Array<Array<double, UserSetup::nEq>, nDeg>>, UserSetup::timeOrder> ti;

  double timestep = 1.0e10;
  Flux<UserSetup::nEq> flux;
  for (int i = mesh->nGhost; i < mesh->Nx - mesh->nGhost; i++) {
    for (int j = mesh->nGhost; j < mesh->Ny - mesh->nGhost; j++) {
      for (int k = mesh->nGhost; k < mesh->Nz - mesh->nGhost; k++) {
        Array<double, UserSetup::nEq> u =
          state.U(mesh_state[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i],
                  0.0, 0.0, 0.0);
        timestep = std::min(timestep, mesh->dx/flux.max_wave_speed_x(u));
        timestep = std::min(timestep, mesh->dy/flux.max_wave_speed_y(u));
        timestep = std::min(timestep, mesh->dz/flux.max_wave_speed_z(u));
      }
    }
  }

  timestep = cf->GetParameter<double>("courant_number")*timestep;

  std::cout << "Time step: " << timestep << "\n";

  std::function<DynArray<Array<Array<double, UserSetup::nEq>, nDeg>>(double, DynArray<Array<Array<double, UserSetup::nEq>, nDeg>>)> L =
                [&rhs](double t,
                       DynArray<Array<Array<double, UserSetup::nEq>, nDeg>> U) -> DynArray<Array<Array<double, UserSetup::nEq>, nDeg>>
    {
    return rhs.Calculate(t, U);
  };

  ti.TakeStep(0.0, timestep, mesh_state, L);

  //double u = 1.0;
  //for (int i = 0; i < 10; i++)
  //  ti.TakeStep(0.0, 0.1, u, [](double t, double u){ return u; });

  //std::cout << u << std::endl;


  delete cf;
  delete mesh;

}

Simulation::~Simulation()
{
  //delete[] state;
}

}
