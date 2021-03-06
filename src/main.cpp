#ifdef USE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <map>
#include <string>

#include "simulation/simulation.hpp"

int main(int argc, char** argv)
{
  // Parse command line arguments
  int nSwitches = 0;                    // Number of command line switches
  double maxWallClockHours = 1.0e10;    // Maximum wallclock hours to run

  // Walk through all command line arguments
  for (int i = 1; i < argc; ++i) {
    std::string arg (argv[i]);

    // Max wall clock hours
    if (arg.compare(std::string("--wallclocklimit")) == 0 ||
        arg.compare(std::string("-wcl")) == 0) {
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

#ifdef USE_MPI
  // Initialize MPI
  MPI_Init(&argc, &argv);
#endif

  // Last argument should be input file name
  char *fileName = argv[argc-1];

  DGHydro::Simulation *s;
  try {
    s = new DGHydro::Simulation(fileName);
  }
  catch (std::exception& e) {
    std::cout << e.what() << '\n';
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return 1;
  }

  delete s;

#ifdef USE_MPI
  MPI_Finalize();
#endif

  return 0;
}
