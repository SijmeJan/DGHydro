#include "mesh.hpp"
#include "../config_file/config_file.hpp"

namespace DGHydro {

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Constructor
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mesh::Mesh(ConfigFile *cf)
{
  nGhost = 1;

  try {
    Nx = cf->GetParameter<int>("Nx") + 2*nGhost;
    startX = nGhost;
    endX = Nx - nGhost;
  }
  catch (std::exception& e) {
    std::cout << e.what() << '\n';
    throw std::runtime_error("No value for Nx, can not create mesh");
  }

  try {
    Ny = cf->GetParameter<int>("Ny") + 2*nGhost;
    startY = nGhost;
    endY = Ny - nGhost;
  }
  catch (std::exception& e) {
    Ny = 1;
    startY = 0;
    endY = 1;
  }

  try {
    Nz = cf->GetParameter<int>("Nz") + 2*nGhost;
    startZ = nGhost;
    endZ = Nz - nGhost;
  }
  catch (std::exception& e) {
    Nz = 1;
    startZ = 0;
    endZ = 1;
  }

  x = new double[Nx];
  y = new double[Ny];
  z = new double[Nz];

  try {
    minX = cf->GetParameter<double>("minX");
    maxX = cf->GetParameter<double>("maxX");

    // Init x: minX and maxX are grid edges
    dx = (maxX - minX)/(double)(Nx - 2*nGhost);
    for (int i = 0; i < Nx; i++)
      x[i] = dx*(double)i + minX - ((double)nGhost - 0.5)*dx;
  }
  catch (std::exception& e) {
    std::cout << e.what() << '\n';
    throw std::runtime_error("No range for x domain, can not create mesh");
  }

  if (Ny > 1) {
    try {
      minY = cf->GetParameter<double>("minY");
      maxY = cf->GetParameter<double>("maxY");

      // Init y: minY and maxY are grid edges
      dy = (maxY - minY)/(double)(Ny - 2*nGhost);
      for (int j = 0; j < Ny; j++)
        y[j] = dy*(double)j + minY - ((double)nGhost - 0.5)*dy;
    }
    catch (std::exception& e) {
      std::cout << e.what() << '\n';
      throw std::runtime_error("No range for y domain, can not create mesh");
    }
  }

  if (Nz > 1) {
    try {
      minZ = cf->GetParameter<double>("minZ");
      maxZ = cf->GetParameter<double>("maxZ");

      // Init z: minZ and maxZ are grid edges
      dz = (maxZ - minZ)/(double)(Nz - 2*nGhost);
      for (int k = 0; k < Nz; k++)
        z[k] = dz*(double)k + minZ - ((double)nGhost - 0.5)*dz;
    }
    catch (std::exception& e) {
      std::cout << e.what() << '\n';
      throw std::runtime_error("No range for z domain, can not create mesh");
    }
  }
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Destructor
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mesh::~Mesh(void)
{
  // Release memory
  delete[] x;
  delete[] y;
  delete[] z;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Domain decomposition
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void Mesh::Decompose(int rank, int num_proc)
{
  if ((Nx - 2*nGhost) % num_proc != 0 ||
      ((Ny - 2*nGhost) % num_proc != 0 && Ny != 1))
    throw std::runtime_error("Numer of cells in x or y not divisible by number of MPI processes");

  // Number of rows residing on this CPU
  int local_nx = (Nx - 2*nGhost)/num_proc + 2*nGhost;

  // Portion of x coordinate local on this process
  double *local_x = new double[local_nx];

  // Init x: xin and xout are grid edges
  for (int i = 0; i < local_nx; i++)
    local_x[i] = x[i] + (double)rank*(double)(local_nx - 2*nGhost)*dx;

  delete[] x;
  x = local_x;
  Nx = local_nx;


}


} // namespace DGHydro
