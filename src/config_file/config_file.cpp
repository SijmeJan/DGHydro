#ifdef USE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include "config_file.hpp"

namespace DGHydro {

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ConfigFile::ConfigFile(char *fileName, int rank)
{
  int warning = 0;

  if (rank == 0) {
    std::ifstream inFile;

    std::cout << "Opening configuration file " << fileName << std::endl;

    // Open configuration file
    inFile.open(fileName);
    if (!inFile.is_open()) {
      warning = 1;
    } else {
      std::string line;
      while (getline(inFile, line)) {
        std::string firstWord, secondWord;

        // Extract first two words from line
        std::istringstream iss(line);
        iss >> firstWord;
        iss >> secondWord;

        // Ignore comments
        if (firstWord.at(0) != '#')
          parameterMap.insert(std::make_pair(firstWord, std::stod(secondWord)));
      }
    }
  }

#ifdef USE_MPI
  // Broadcast warning to see if everything is OK
  MPI_Bcast(&warning, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  if (warning != 0)
    throw std::runtime_error("Error opening configuration file");

  // Broadcast size of map to all processes
  int mapSize = parameterMap.size();
#ifdef USE_MPI
  MPI_Bcast(&mapSize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  // Iterate through all elements
  std::map<std::string, double>::iterator it = parameterMap.begin();
  for (int i = 0; i < mapSize; i++) {
    std::string line;
    int line_size = 0;
    if (rank == 0) {
      line = it->first;
      line_size = line.size();
    }
#ifdef USE_MPI
    MPI_Bcast(&line_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    if (rank != 0)
      line.resize(line_size);

#ifdef USE_MPI
    MPI_Bcast(const_cast<char*>(line.data()), line_size,
              MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

    double value = 0.0;
    if (rank == 0)
      value = it->second;

#ifdef USE_MPI
    MPI_Bcast(&value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    if (rank != 0)
      parameterMap.insert(std::make_pair(line, value));

    if (rank == 0)
      it++;
  }
}

ConfigFile::~ConfigFile(void)
{
  // Nothing to be done
}

void ConfigFile::List()
{
  // Iterate through all elements
  std::map<std::string, double>::iterator it = parameterMap.begin();
  while(it != parameterMap.end()) {
    std::cout << it->first << " :: " << it->second << std::endl;
    it++;
  }
}

}
