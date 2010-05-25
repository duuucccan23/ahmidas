#include "Propagator.ih"

#include <L1/Smear/Jacobi.h>

namespace Core
{
  Propagator &Propagator::smearJacobi(double const kappa, size_t const iterations, Field< QCD::Gauge > &gauge_field)
  {
    isolate();
    Smear::Jacobi JacobiTool(kappa);
    JacobiTool.smear(d_components, gauge_field, iterations);
    return *this;
  }

  Propagator &Propagator::smearJacobi(double const kappa, size_t const iterations, Field< QCD::Gauge > &gauge_field, size_t const timeslice)
  {
    isolate();
    Smear::Jacobi JacobiTool(kappa);
    JacobiTool.smear(d_components, gauge_field, iterations, timeslice);
    return *this;
  }
}
