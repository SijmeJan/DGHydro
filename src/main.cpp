#include <mpi.h>

#include <iostream>
#include <map>
#include <string>

//#include "mesh/mesh.hpp"
//#include "config_file/config_file.hpp"
#include "simulation/simulation.hpp"

int main(int argc, char** argv)
{
  // Parse command line arguments
  int nSwitches = 0;                    // Number of command line switches
  double maxWallClockHours = 1.0e10;    // Maximum wallclock hours to run

  // Walk through all command line arguments
  for (int i = 1; i < argc; ++i) {
    // Max wall clock hours
    if (strcmp(argv[i], "--wallclocklimit") == 0 ||
        strcmp(argv[i], "-wcl") == 0) {
      maxWallClockHours = atof(argv[i+1]);
      std::cout << "Maximum wall clock time: " << maxWallClockHours
                << " hours" << std::endl;
      nSwitches += 2;
    }
  }

  // Check for correct number of arguments
  if (argc != 2 + nSwitches) {
    std::cout << "Usage: " << argv[0]
              << " [-wcl maxWallClockHours]"
              << " filename"
              << std::endl;
    return 1;
  }

  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Last argument should be input file name
  char *fileName = argv[argc-1];

  try {
    DGHydro::Simulation *s = new DGHydro::Simulation(fileName);
  }
  catch (std::exception& e) {
    std::cout << e.what() << '\n';
    MPI_Finalize();
    return 1;
  }

  MPI_Finalize();

  return 0;
}
