#ifdef USE_MPI
#include <mpi.h>
#endif
#include <iostream>
#include <fstream>

#include "../array/array.hpp"
#include "../array/dynarray.hpp"
#include "../state/state.hpp"

#include "simulation.hpp"
#include "../config_file/config_file.hpp"
#include "../mesh/mesh.hpp"
#include "../initial/initial.hpp"
#include "../rhs/rhs.hpp"
#include "../flux/flux.hpp"
#include "time_integrator.hpp"

namespace DGHydro {

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Constructor
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Simulation::Simulation(char *fileName)
{
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

  if (rank == 0)
    std::cout << "Number of space dimensions: " << UserSetup::nDim
              << "\nNumber of equations: " << UserSetup::nEq
              << "\nSpace order: " << UserSetup::maxOrder + 1
              << "\nTime order: " << UserSetup::timeOrder
              << std::endl;

  // Read configuration file
  try {
    cf = new ConfigFile(fileName, rank);
    cf->List();
  }
  catch (...) {
    std::cout << "Error reading configuration file!" << std::endl;
    throw;
  }

  // Build mesh
  try {
    mesh = new Mesh(cf);
    mesh->Decompose(rank, num_proc);
  }
  catch (std::exception& e) {
    std::cout << e.what() << '\n';
    throw std::runtime_error("Could not create simulation");
  }

  if (mesh->Ny == 1 && UserSetup::nDim > 1) {
    delete mesh;
    delete cf;
    throw std::runtime_error("Need Ny in config file for multidimensional simulation");
  }
  if (mesh->Nz == 1 && UserSetup::nDim == 3) {
    delete mesh;
    delete cf;
    throw std::runtime_error("Need Nz in config file for 3-dimensional simulation");
  }

  state = new State<UserSetup::nEq,
                    UserSetup::maxOrder,
                    UserSetup::nDim>(mesh);

  mesh_state = new DynArray<t_state_deg>(mesh->Nx*mesh->Ny*mesh->Nz);
  mesh_state[0] = 0.0;

  cfl = cf->GetParameter<double>("courant_number");

  std:: cout << "Restoring from dump...\n";

  try {
    RestoreFromDump(0);
  }
  catch (std::exception& e) {
    std::cout << e.what() << '\n';
    throw std::runtime_error("Could not create simulation");
  }

  std:: cout << "Initial conditions...\n";

  // Set initial conditions
  InitialConditions<UserSetup::nEq> ic(cf);
  for (int i = mesh->startX; i < mesh->endX; i++)
    for (int j = mesh->startY; j < mesh->endY; j++)
      for (int k = mesh->startZ; k < mesh->endZ; k++)
        mesh_state[0][k*mesh->Nx*mesh->Ny + j*mesh->Nx + i] =
          state->DoF(i, j, k, ic);

  std:: cout << "Saving to dump...\n";

  try {
    Dump(0);
  }
  catch (std::exception& e) {
    std::cout << e.what() << '\n';
    throw std::runtime_error("Could not create simulation");
  }

  std::cout << "Done\n";

  /*
  try {
    RestoreFromDump(0);
  }
  catch (std::exception& e) {
    std::cout << e.what() << '\n';
    throw std::runtime_error("Could not create simulation");
  }
  */
  std::cout << "Watte\n";

  // Set up right-hand side of d_t U = RHS(t, U)
  RightHandSide<UserSetup::nEq, UserSetup::maxOrder, UserSetup::nDim> rhs(mesh);

  std::cout << "Hallo" << std::endl;

  // Set up time integrator
  TimeIntegrator<DynArray<t_state_deg>, UserSetup::timeOrder> ti;

  std::cout << "Doei\n";

  std::function<DynArray<t_state_deg>(double, DynArray<t_state_deg>)>
    L = [&rhs](double t, DynArray<t_state_deg> U) -> DynArray<t_state_deg> {
    return rhs.Calculate(t, U);
  };

  std::cout << "Starting time loop..." << std::endl;

  double time = 0.0;
  while (time < 1.0) {
    double timestep = CalcTimeStep();

    ti.TakeStep(time, timestep, mesh_state[0], L);

    time += timestep;

    std::cout << time << std::endl;
  }

}

Simulation::~Simulation()
{
  delete state;
  delete mesh_state;
  delete cf;
  delete mesh;
}

double Simulation::CalcTimeStep()
{
  double timestep = 1.0e10;
  Flux<UserSetup::nEq> flux;
  for (int i = mesh->startX; i < mesh->endX; i++) {
    for (int j = mesh->startY; j < mesh->endY; j++) {
      for (int k = mesh->startZ; k < mesh->endZ; k++) {
        std::cout << "Times step " << i << " " << j << " " << k << "\n";
        t_state u =
          state->U(mesh_state[0][k*mesh->Nx*mesh->Ny + j*mesh->Nx + i],
                   0.0, 0.0, 0.0);
        timestep = std::min(timestep, mesh->dx/flux.max_wave_speed_x(u));
        timestep = std::min(timestep, mesh->dy/flux.max_wave_speed_y(u));
        timestep = std::min(timestep, mesh->dz/flux.max_wave_speed_z(u));
      }
    }
  }

  return cfl*timestep;
}

void Simulation::Dump(int nDump)
{
  char fname[15];
  sprintf(fname, "dump%3.3d.dat", nDump);

  std::ofstream wf(fname, std::ios::out | std::ios::binary);
  if (!wf)
    throw std::runtime_error("Can not open output file for writing!");

  std::cout << "Size of state: " << sizeof(t_state_deg)
            << " " << sizeof(double)
            << " " << UserSetup::nDeg << std::endl;
  for (int i = mesh->startX; i < mesh->endX; i++)
    for (int j = mesh->startY; j < mesh->endY; j++)
      for (int k = mesh->startZ; k < mesh->endZ; k++)
        for (int n = 0; n < UserSetup::nDeg; n++)
          wf.write((char *) &mesh_state[0][k*mesh->Nx*mesh->Ny + j*mesh->Nx + i][n],
                 sizeof(t_state));

  wf.close();
  if (!wf.good())
    throw std::runtime_error("Could not write output file!");
}

void Simulation::RestoreFromDump(int nDump)
{
  char fname[15];
  sprintf(fname, "dump%3.3d.dat", nDump);

  std::ifstream rf(fname, std::ios::out | std::ios::binary);
  if (!rf)
    throw std::runtime_error("Can not open dump file for reading!");

  int index = 0;
  for (int i = mesh->startX; i < mesh->endX; i++)
    for (int j = mesh->startY; j < mesh->endY; j++)
      for (int k = mesh->startZ; k < mesh->endZ; k++)
        for (int n = 0; n < UserSetup::nDeg; n++) {
          rf.read((char *) &mesh_state[0][k*mesh->Nx*mesh->Ny + j*mesh->Nx + i][n],
                  sizeof(t_state));
          std::cout << "Read index " << index++ << std::endl;
        }

  rf.close();
  if (!rf.good())
    throw std::runtime_error("Could not read dump input file!");
}

}
