#include "Propagator.ih"

namespace Core
{
  template< >
  StochasticPropagator< 4 >::StochasticPropagator(StochasticPropagator< 4 > const &other, size_t const timeslice)
  : Propagator(dynamic_cast< Propagator const& >(other), timeslice)
  {}
}

