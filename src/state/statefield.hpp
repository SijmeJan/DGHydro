#ifndef DG_STATEFIELD_HPP
#define DG_STATEFIELD_HPP


namespace DGHydro {

  template<int nEq, int nDeg>
  class StateField {
  public:
    StateField(int _Nx, int _Ny, int _Nz) {
      Nx = _Nx;
      Ny = _Ny;
      Nz = _Nz;

      data = new double[nDeg*nEq*Nx*Ny*Nz];
    };
    ~StateField() {
      delete[] data;
    };

    // Get state at equation number N, d.o.f number M
    template<int N, int M>
    double& get(int i) {
      return data[M*nEq*Nx + N*Nx + i];
    }


  private:
    int Nx, Ny, Nz;
    double *data;
  };

}

#endif // DG_STATEFIELD_HPP
