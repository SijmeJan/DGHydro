#ifndef DG_RHS_HPP
#define DG_RHS_HPP

#include "../array/array.hpp"
#include "../flux/flux.hpp"
#include "../state/basis.hpp"
#include "../mesh/mesh.hpp"
#include "temp.hpp"

namespace DGHydro {

  template<int nEq, int maxOrder, int nDim>
  class RightHandSide {
  public:
    RightHandSide(Mesh *mesh) : mesh(mesh) {
    };
    ~RightHandSide() {
    };

    // Number of degrees of freedom
    const static int nDeg =
      (nDim == 1)*(maxOrder + 1) +
      (nDim == 2)*(maxOrder + 1)*(maxOrder + 2)/2 +
      (nDim == 3)*(maxOrder + 1)*(maxOrder + 2)*(maxOrder + 3)/6;

    using t_state_deg = Array<Array<double, nEq>, nDeg>;
    using t_state = Array<double, nEq>;

    DynArray<t_state_deg>
    Calculate(double t, DynArray<t_state_deg>& U) {
      DynArray<t_state_deg> data =
        DynArray<t_state_deg>(mesh->Nx*mesh->Ny*mesh->Nz);
      data = 0.0;

      SetBoundaries(U);

      for (int i = mesh->startX; i < mesh->endX; i++)
        for (int j = mesh->startY; j < mesh->endY; j++)
          for (int k = mesh->startZ; k < mesh->endZ; k++)
            data[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i] =
              VolumeFluxIntegral(U[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i]);

      for (int i = mesh->startX; i < mesh->endX; i++)
        for (int j = mesh->startY; j < mesh->endY; j++)
          for (int k = mesh->startZ; k < mesh->endZ; k++)
            data[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i] +=
              SurfaceFluxIntegralX(U[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i],
                                   U[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i - 1],
                                   U[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i + 1]);

      if (nDim > 1)
        for (int i = mesh->startX; i < mesh->endX; i++)
          for (int j = mesh->startY; j < mesh->endY; j++)
            for (int k = mesh->startZ; k < mesh->endZ; k++)
              data[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i] +=
                SurfaceFluxIntegralY(U[k*mesh->Nx*mesh->Ny+ j*mesh->Nx + i],
                                     U[k*mesh->Nx*mesh->Ny+(j-1)*mesh->Nx + i],
                                     U[k*mesh->Nx*mesh->Ny+(j+1)*mesh->Nx + i]);

      if (nDim > 2)
        for (int i = mesh->startX; i < mesh->endX; i++)
          for (int j = mesh->startY; j < mesh->endY; j++)
            for (int k = mesh->startZ; k < mesh->endZ; k++)
              data[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i] +=
                SurfaceFluxIntegralZ(U[k*mesh->Nx*mesh->Ny + j*mesh->Nx + i],
                                     U[(k-1)*mesh->Nx*mesh->Ny+j*mesh->Nx + i],
                                     U[(k+1)*mesh->Nx*mesh->Ny+j*mesh->Nx + i]);

      return data;
    }


  private:
    Flux<nEq> flux;
    BasisFunctions<nDim, maxOrder> bf;
    CubeIntegral<t_state, maxOrder + 1> ci;

    Mesh *mesh;

    // Calculate state from degrees of freedom
    t_state U(t_state_deg& s, double x, double y, double z) {
      t_state result;
      result = s[0]*bf(0, x, y, z);

      // Sum of components times basis function
      for (int j = 1; j < nDeg; j++)
        result += s[j]*bf(j, x, y, z);

      result *= (1 << nDim);

      return result;
    }

    // LLF interface flux x
    t_state Fint_x(t_state_deg& sL, t_state_deg& sR,
                   double y, double z) {
      t_state UL = U(sL, 1, y, z);
      t_state UR = U(sR, -1, y, z);

      double lambda = std::max(flux.max_wave_speed_x(UL),
                               flux.max_wave_speed_x(UR));
      return 0.5*(flux.x(UR) + flux.x(UL) - lambda*(UR - UL));
    }
    // LLF interface flux y
    t_state Fint_y(t_state_deg& sL, t_state_deg& sR,
                   double x, double z) {
      t_state UL = U(sL, x, 1, z);
      t_state UR = U(sR, x, -1, z);

      double lambda = std::max(flux.max_wave_speed_y(UL),
                               flux.max_wave_speed_y(UR));
      return 0.5*(flux.y(UR) + flux.y(UL) - lambda*(UR - UL));
    }
    // LLF interface flux z
    t_state Fint_z(t_state_deg& sL, t_state_deg& sR,
                   double x, double y) {
      t_state UL = U(sL, x, y, 1);
      t_state UR = U(sR, x, y, -1);

      double lambda = std::max(flux.max_wave_speed_z(UL),
                               flux.max_wave_speed_z(UR));
      return 0.5*(flux.z(UR) + flux.z(UL) - lambda*(UR - UL));
    }


    t_state_deg SurfaceFluxIntegralX(t_state_deg& s,
                                     t_state_deg& sL,
                                     t_state_deg& sR) {
      t_state_deg result(0.0);

      for (int j = 0; j < nDeg; j++) {
        // 3D: integrate over 2D surface
        if (nDim == 3) {
          // Function to integrate
          std::function<t_state(double, double)> f =
            [&s, &sL, &sR, this, j] (double y, double z) -> t_state {

            return 0.25*(Fint_x(sL, s, y, z)*bf(j, -1, y, z) -
                         Fint_x(s, sR, y, z)*bf(j, 1, y, z))/mesh->dx;
          };

          result[j] = ci.vol2d(f);
        }
        // 2D: integrate over 1D line
        if (nDim == 2) {
          // Function to integrate
          std::function<t_state(double)> f =
            [&s, &sL, &sR, this, j] (double y) -> t_state {

            return 0.5*(Fint_x(sL, s, y, 0)*bf(j, -1, y, 0) -
                         Fint_x(s, sR, y, 0)*bf(j, 1, y, 0))/mesh->dx;
          };

          result[j] = ci.vol1d(f);
        }
        // 1D: no integration necessary
        if (nDim == 1) {
          result[j] = (Fint_x(sL, s, 0, 0)*bf(j, -1, 0, 0) -
                       Fint_x(s, sR, 0, 0)*bf(j, 1, 0, 0))/mesh->dx;
        }
      }

      return result;
    };

    t_state_deg SurfaceFluxIntegralY(t_state_deg& s,
                                     t_state_deg& sL,
                                     t_state_deg& sR) {
      t_state_deg result(0.0);

      for (int j = 0; j < nDeg; j++) {
        if (nDim == 3) {
          // Function to integrate
          std::function<t_state(double, double)> f =
            [&s, &sL, &sR, this, j] (double x, double z) -> t_state {

            return 0.25*(Fint_x(sL, s, x, z)*bf(j, x, -1, z) -
                         Fint_x(s, sR, x, z)*bf(j, x, 1, z))/mesh->dy;
          };

          result[j] = ci.vol2d(f);
        }
        if (nDim == 2) {
          // Function to integrate
          std::function<t_state(double)> f =
            [&s, &sL, &sR, this, j] (double x) -> t_state {

            return 0.5*(Fint_x(sL, s, x, 0)*bf(j, x, -1, 0) -
                        Fint_x(s, sR, x, 0)*bf(j, x, 1, 0))/mesh->dy;
          };

          result[j] = ci.vol1d(f);
        }
      }

      return result;
    };
    t_state_deg SurfaceFluxIntegralZ(t_state_deg& s,
                                     t_state_deg& sL,
                                     t_state_deg& sR) {
      t_state_deg result(0.0);

      for (int j = 0; j < nDeg; j++) {
        // Function to integrate
        std::function<t_state(double, double)> f =
          [&s, &sL, &sR, this, j] (double x, double y) -> t_state {

          return 0.25*(Fint_x(sL, s, x, y)*bf(j, x, y, -1) -
                       Fint_x(s, sR, x, y)*bf(j, x, y, 1))/mesh->dz;
        };

        result[j] = ci.vol2d(f);
      }

      return result;
    };

    t_state_deg VolumeFluxIntegral(t_state_deg& s) {
      t_state_deg result(0.0);

      for (int j = 0; j < nDeg; j++) {
        if (nDim == 3) {
          // Function to integrate
          std::function<t_state(double, double, double)> f =
            [&s, this, j] (double x, double y, double z) -> t_state {

            double a = 0.25*bf.x_derivative(j, x, y, z)/mesh->dx;
            double b = 0.25*bf.y_derivative(j, x, y, z)/mesh->dy;
            double c = 0.25*bf.z_derivative(j, x, y, z)/mesh->dz;

            t_state u = U(s, x, y, 0);
            return flux.x(u)*a + flux.y(u)*b + flux.z(u)*c;
          };

          result[j] = ci.vol3d(f);
        }
        if (nDim == 2) {
          // Function to integrate
          std::function<t_state(double, double)> f =
            [&s, this, j] (double x, double y) -> t_state {

            double a = 0.5*bf.x_derivative(j, x, y, 0)/mesh->dx;
            double b = 0.5*bf.y_derivative(j, x, y, 0)/mesh->dy;

            t_state u = U(s, x, y, 0);
            return flux.x(u)*a + flux.y(u)*b;
          };

          result[j] = ci.vol2d(f);
        }
        if (nDim == 1) {
          // Function to integrate
          std::function<t_state(double)> f =
            [&s, this, j] (double x) -> t_state {

            double a = bf.x_derivative(j, x, 0, 0)/mesh->dx;

            t_state u = U(s, x, 0, 0);
            return flux.x(u)*a;
          };

          result[j] = ci.vol1d(f);
        }
      }

      return result;
    };

    void SetBoundaries(DynArray<t_state_deg>& U) {
      // Periodic x boundary
      for (int j = mesh->startY; j < mesh->endY; j++) {
        for (int k = mesh->startZ; k < mesh->endZ; k++) {
          U[k*mesh->Nx*mesh->Ny + j*mesh->Nx + 0] =
            U[k*mesh->Nx*mesh->Ny + j*mesh->Nx + mesh->Nx - 2];

          U[k*mesh->Nx*mesh->Ny + j*mesh->Nx + mesh->Nx - 1] =
            U[k*mesh->Nx*mesh->Ny + j*mesh->Nx + 1];
        }
      }

      if (nDim > 1) {
        for (int i = mesh->startX; i < mesh->endX; i++) {
          for (int k = mesh->startZ; k < mesh->endZ; k++) {
            U[k*mesh->Nx*mesh->Ny + 0*mesh->Nx + i] =
              U[k*mesh->Nx*mesh->Ny + (mesh->Ny-2)*mesh->Nx + i];

            U[k*mesh->Nx*mesh->Ny + (mesh->Ny-1)*mesh->Nx + i] =
              U[k*mesh->Nx*mesh->Ny + 1*mesh->Nx + i];
          }
        }
      }

      if (nDim > 2) {
        for (int i = mesh->startX; i < mesh->endX; i++) {
          for (int j = mesh->startY; j < mesh->endY; j++) {
            U[0*mesh->Nx*mesh->Ny + j*mesh->Nx + i] =
              U[(mesh->Nz-2)*mesh->Nx*mesh->Ny + j*mesh->Nx + i];

            U[(mesh->Nz-1)*mesh->Nx*mesh->Ny + j*mesh->Nx + i] =
              U[1*mesh->Nx*mesh->Ny + j*mesh->Nx + i];
          }
        }
      }
    }

  };

} // namespace DGHydro

#endif  // DG_RHS_HPP
