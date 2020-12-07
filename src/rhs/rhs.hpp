#ifndef DG_RHS_HPP
#define DG_RHS_HPP

#include "../state/array.hpp"
#include "../state/mesharray.hpp"
#include "../flux/flux.hpp"
#include "../state/basis.hpp"
#include "../mesh/mesh.hpp"
#include "temp.hpp"

namespace DGHydro {

  template<int nEq, int maxOrder, int nDim>
  class RightHandSide {
  public:
    RightHandSide(Mesh *mesh) : mesh(mesh) {
      //data = new Array<Array<double, nEq>, nDeg>[mesh->Nx*mesh->Ny*mesh->Nz];
      data = new MeshArray<Array<Array<double, nEq>, nDeg>>(mesh->Nx,
                                                            mesh->Ny,
                                                            mesh->Nz);
    };
    ~RightHandSide() {
      delete data;
    };

    // Number of degrees of freedom
    const static int nDeg =
      (nDim == 1)*maxOrder +
      (nDim == 2)*(maxOrder + 1)*(maxOrder + 2)/2 +
      (nDim == 3)*(maxOrder + 1)*(maxOrder + 2)*(maxOrder + 3)/6;

    //void Calculate(Array<Array<double, nEq>, nDeg> *U) {
    void Calculate(MeshArray<Array<Array<double, nEq>, nDeg>>& U) {
      /*
      for (int i = mesh->nGhost; i < mesh->Nx - mesh->nGhost; i++)
        for (int j = mesh->nGhost; j < mesh->Ny - mesh->nGhost; j++)
          for (int k = mesh->nGhost; k < mesh->Nz - mesh->nGhost; k++)
            data[0][k*mesh->Nx*mesh->Ny + j*mesh->Nx + i] =
              VolumeFluxIntegral(U[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i]);

      for (int i = mesh->nGhost; i < mesh->Nx - mesh->nGhost; i++) {
        for (int j = mesh->nGhost; j < mesh->Ny - mesh->nGhost; j++) {
          for (int k = mesh->nGhost; k < mesh->Nz - mesh->nGhost; k++) {
            data[0][k*mesh->Nx*mesh->Ny + j*mesh->Nx + i] +=
              SurfaceFluxIntegralX(U[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i],
                                   U[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i - 1],
                                   U[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i + 1]);

            data[0][k*mesh->Nx*mesh->Ny + j*mesh->Nx + i] +=
              SurfaceFluxIntegralY(U[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i],
                                   U[k*mesh->Nx*mesh->Ny + (j - 1)*mesh->Nx + i],
                                   U[k*mesh->Nx*mesh->Ny + (j + 1)*mesh->Nx + i]);

            data[0][k*mesh->Nx*mesh->Ny + j*mesh->Nx + i] +=
              SurfaceFluxIntegralZ(U[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i],
                                   U[(k - 1)*mesh->Nx*mesh->Ny + j*mesh->Nx + i],
                                   U[(k + 1)*mesh->Nx*mesh->Ny + j*mesh->Nx + i]);
          }
        }
      }
      */
    }

    Array<Array<double, nEq>, nDeg> SurfaceFluxIntegralX(Array<Array<double, nEq>, nDeg>& s,
                                                         Array<Array<double, nEq>, nDeg>& sL,
                                                         Array<Array<double, nEq>, nDeg>& sR) {
      Array<Array<double, nEq>, nDeg> result;

      for (int j = 0; j < nDeg; j++) {
        // Function to integrate
        std::function<Array<double, nEq>(double, double)> f =
          [&s, &sL, &sR, this, j] (double y, double z) -> Array<double, nEq> {

          return 0.25*(Fint_x(sL, s, y, z)*bf(j, -1, y, z) -
                       Fint_x(s, sR, y, z)*bf(j, 1, y, z))/mesh->dx;
        };

        result[j] = ci.vol2d(f);
      }

      return result;
    };
    Array<Array<double, nEq>, nDeg> SurfaceFluxIntegralY(Array<Array<double, nEq>, nDeg>& s,
                                                         Array<Array<double, nEq>, nDeg>& sL,
                                                         Array<Array<double, nEq>, nDeg>& sR) {
      Array<Array<double, nEq>, nDeg> result;

      for (int j = 0; j < nDeg; j++) {
        // Function to integrate
        std::function<Array<double, nEq>(double, double)> f =
          [&s, &sL, &sR, this, j] (double x, double z) -> Array<double, nEq> {

          return 0.25*(Fint_x(sL, s, x, z)*bf(j, x, -1, z) -
                       Fint_x(s, sR, x, z)*bf(j, x, 1, z))/mesh->dy;
        };

        result[j] = ci.vol2d(f);
      }

      return result;
    };
    Array<Array<double, nEq>, nDeg> SurfaceFluxIntegralZ(Array<Array<double, nEq>, nDeg>& s,
                                                         Array<Array<double, nEq>, nDeg>& sL,
                                                         Array<Array<double, nEq>, nDeg>& sR) {
      Array<Array<double, nEq>, nDeg> result;

      for (int j = 0; j < nDeg; j++) {
        // Function to integrate
        std::function<Array<double, nEq>(double, double)> f =
          [&s, &sL, &sR, this, j] (double x, double y) -> Array<double, nEq> {

          return 0.25*(Fint_x(sL, s, x, y)*bf(j, x, y, -1) -
                       Fint_x(s, sR, x, y)*bf(j, x, y, 1))/mesh->dz;
        };

        result[j] = ci.vol2d(f);
      }

      return result;
    };

    Array<Array<double, nEq>, nDeg> VolumeFluxIntegral(Array<Array<double, nEq>, nDeg>& s) {
      Array<Array<double, nEq>, nDeg> result;

      for (int j = 0; j < nDeg; j++) {
        // Function to integrate
        std::function<Array<double, nEq>(double, double, double)> f =
          [&s, this, j] (double x, double y, double z) -> Array<double, nEq> {

          return 0.25*flux.x(U(s, x, y, z))*bf.x_derivative(j, x, y, z)/mesh->dx +
                 0.25*flux.y(U(s, x, y, z))*bf.y_derivative(j, x, y, z)/mesh->dy +
                 0.25*flux.z(U(s, x, y, z))*bf.z_derivative(j, x, y, z)/mesh->dz;
        };

        result[j] = ci.vol3d(f);
      }

      return result;
    };

    //Array<Array<double, nEq>, nDeg> *data;
    MeshArray<Array<Array<double, UserSetup::nEq>, nDeg>> *data;
  private:
    Flux<nEq> flux;
    BasisFunctions<nDim, maxOrder> bf;
    CubeIntegral<Array<double, nEq>, maxOrder + 1> ci;

    Mesh *mesh;

    // Calculate state from degrees of freedom
    Array<double, nEq> U(Array<Array<double, nEq>, nDeg>& s, double x, double y, double z) {
      Array<double, nEq> result;
      result = 0.0;

      // Sum of components times basis function
      for (int j = 0; j < nDeg; j++)
        result += s[j]*bf(j, x, y, z);

      return result*(1 << nDim);
    }

    // LLF interface flux x
    Array<double, nEq> Fint_x(Array<Array<double, nEq>, nDeg>& sL,
                              Array<Array<double, nEq>, nDeg>& sR,
                              double y, double z) {
      Array<double, nEq> UL = U(sL, 1, y, z);
      Array<double, nEq> UR = U(sR, -1, y, z);

      double lambda = std::max(flux.max_wave_speed_x(UL),
                               flux.max_wave_speed_x(UR));
      return 0.5*(flux.x(UR) + flux.x(UL) - lambda*(UR - UL));
    }
    // LLF interface flux y
    Array<double, nEq> Fint_y(Array<Array<double, nEq>, nDeg>& sL,
                              Array<Array<double, nEq>, nDeg>& sR,
                              double x, double z) {
      Array<double, nEq> UL = U(sL, x, 1, z);
      Array<double, nEq> UR = U(sR, x, -1, z);

      double lambda = std::max(flux.max_wave_speed_y(UL),
                               flux.max_wave_speed_y(UR));
      return 0.5*(flux.y(UR) + flux.y(UL) - lambda*(UR - UL));
    }
    // LLF interface flux z
    Array<double, nEq> Fint_z(Array<Array<double, nEq>, nDeg>& sL,
                              Array<Array<double, nEq>, nDeg>& sR,
                              double x, double y) {
      Array<double, nEq> UL = U(sL, x, y, 1);
      Array<double, nEq> UR = U(sR, x, y, -1);

      double lambda = std::max(flux.max_wave_speed_z(UL),
                               flux.max_wave_speed_z(UR));
      return 0.5*(flux.z(UR) + flux.z(UL) - lambda*(UR - UL));
    }

  };

} // namespace DGHydro

#endif  // DG_RHS_HPP
