#ifndef DG_TIME_HPP
#define DG_TIME_HPP

namespace DGHydro {
  template<class T, int order>
  class TimeIntegrator {
  public:
    TimeIntegrator() {
    };
    ~TimeIntegrator() {
    };

    void TakeStep(double start_time,
                  double timestep, T& U,
                  std::function<T(double, T)> L)
    {
      if (order == 1)
        U += L(start_time, U)*timestep;

      if (order == 2) {
        T U1 = U + timestep*L(start_time, U);
        U = 0.5*U + 0.5*U1 + 0.5*timestep*L(start_time + 0.5*timestep, U1);
      }

      if (order == 3) {
        T U1 = U + timestep*L(start_time, U);
        U1 = 0.75*U + 0.25*U1 + 0.25*timestep*L(start_time + timestep, U1);
        U = U/3 + 2*U1/3 + 2*timestep*L(start_time + 0.5*timestep, U1)/3;
      }

    };

  private:
  };

} // namespace DGHydro

#endif  // DG_TIME_HPP
