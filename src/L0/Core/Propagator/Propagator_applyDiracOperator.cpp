#include "Propagator.ih"

namespace Core
{
  // applies the (twisted Mass) Dirac operator with fixed boundary conditions
  // (temporal gauge links starting at timeslice t_boundary pick up a factor -1)
  Propagator Propagator::applyDiracOperator(Field< QCD::Gauge > const &gauge_field, double const kappa, size_t const t_boundary)
  {
    Propagator result(*this);
    Propagator neighbours(*this);

    // nothing implemented so far :-(


    return result;
  }
}