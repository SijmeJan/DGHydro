#ifndef DG_STATE_HPP
#define DG_STATE_HPP

#include "../mesh/mesh.hpp"
#include "../state/basis.hpp"
#include "../gauss_quad/integral.hpp"

namespace DGHydro {

  template<int nEq, int maxOrder, int nDim>
  class State {
  public:
    State(Mesh *mesh) : mesh(mesh) {
    };
    ~State() {
    };

    // Number of degrees of freedom
    const static int nDeg =
      (nDim == 1)*(maxOrder + 1) +
      (nDim == 2)*(maxOrder + 1)*(maxOrder + 2)/2 +
      (nDim == 3)*(maxOrder + 1)*(maxOrder + 2)*(maxOrder + 3)/6;

    BasisFunctions<nDim, maxOrder> bf;
    CubeIntegral<Array<double, nEq>, maxOrder + 1> ci;

    // Calculate state from degrees of freedom
    Array<double, nEq> U(Array<Array<double, nEq>, nDeg>& s,
                         double x, double y, double z) {
      Array<double, nEq> result;
      result = s[0]*bf(0, x, y, z);

      // Sum of components times basis function
      for (int j = 1; j < nDeg; j++)
        result += s[j]*bf(j, x, y, z);

      result *= (1 << nDim);

      return result;
    }

    // Calculate degrees of freedom from state
    Array<Array<double, nEq>, nDeg>
    DoF(int i, int j, int k,
        std::function<Array<double, nEq>(double, double, double)> func) {
      Array<Array<double, nEq>, nDeg> result(0.0);

      for (int n = 0; n < nDeg; n++) {
        // Function to integrate
        std::function<Array<double, nEq>(double, double, double)> f =
          [this, func, n, i, j, k] (double x, double y, double z) -> Array<double, nEq> {

          return func(mesh->x[i] + 0.5*mesh->dx*x,
                      mesh->y[j] + 0.5*mesh->dy*y,
                      mesh->z[k] + 0.5*mesh->dz*z)*bf(n, x, y, z);
        };

        result[n] = ci.vol3d(f);
      }

      result /= (1 << nDim);

      return result;
    }

  private:
    Mesh *mesh;
  };

}

#endif // DG_STATE_HPP
