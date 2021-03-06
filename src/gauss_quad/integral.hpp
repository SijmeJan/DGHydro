#ifndef DG_INTEGRAL_HPP
#define DG_INTEGRAL_HPP

#include "gauss_quad.hpp"

namespace DGHydro {

  template <class T, int nPoint>
  class CubeIntegral {
  public:
    CubeIntegral() {
      gq = new GaussQuad(nPoint);
    }
    ~CubeIntegral(void) {
      delete gq;
    }

    T vol3d(std::function<T(double, double, double)> func) {
      T res(0.0);         // Allocate 1 T

      for (int k = 0; k < nPoint; k++)
        for (int j = 0; j < nPoint; j++)
          for (int i = 0; i < nPoint; i++)
            res += (gq->w[k]*gq->w[j]*gq->w[i])*
              func(gq->x[i], gq->x[j], gq->x[k]);

      return res;
    };

    T vol2d(std::function<T(double, double)> func) {
      T res(0.0);

      for (int j = 0; j < nPoint; j++)
        for (int i = 0; i < nPoint; i++)
          res += (gq->w[j]*gq->w[i])*func(gq->x[i], gq->x[j]);

      return res;
    };

    T vol1d(std::function<T(double)> func) {
      T res(0.0);

      for (int i = 0; i < nPoint; i++)
        res += gq->w[i]*func(gq->x[i]);

      return res;
    };

  private:
    GaussQuad *gq;
  };


} // namespace DGHydro

#endif  // DG_INTEGRAL_HPP
