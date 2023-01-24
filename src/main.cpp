//#include <Kokkos_Core.hpp>

int main( int argc, char* argv[] ) {
  bool initKokkosBeforeMPI = false;

  // When running on GPUS with Omnipath network,
  // Kokkos needs to be initialised *before* the MPI layer
#ifdef KOKKOS_ENABLE_CUDA
  if(std::getenv("PSM2_CUDA") != NULL) {
    initKokkosBeforeMPI = true;
  }
#endif

  //if(initKokkosBeforeMPI)  Kokkos::initialize( argc, argv );

#ifdef WITH_MPI
  MPI_Init(&argc,&argv);
#endif

  //if(!initKokkosBeforeMPI) Kokkos::initialize( argc, argv );

  return 0;
}
