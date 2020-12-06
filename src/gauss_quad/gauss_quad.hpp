#ifndef DG_GAUSS_HPP
#define DG_GAUSS_HPP

namespace DGHydro {

  class GaussQuad {
  public:
    GaussQuad(int _n);
    ~GaussQuad(void);

    int n;
    double *x, *w;

  private:
    void FindAbscissae();
    void FindWeights();

    template <typename T> int sgn(T val) {
      return (T(0) < val) - (val < T(0));
    }
  };

} // namespace DGHydro

#endif  // DG_GAUSS_HPP
