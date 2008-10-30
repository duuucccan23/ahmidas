#ifndef GUARD_SOURCE_STOCHASTIC_H
#define GUARD_SOURCE_STOCHASTIC_H

#include <L0/Core/Field.h>
#include <L0/SU3/Vector.h>
#include <L0/Base/Random.h>
#include <cmath>

namespace Source
{
  template< size_t L, size_t T >
  class Stochastic
  {
    size_t const                      d_timeSlice;
    Core::Field< SU3::Vector, L, 1 >  d_source;

    public:
      Stochastic();
      size_t timeSlice() const;
      double const *source() const;
  };
}

#include "Stochastic/Stochastic.inlines"
#include "Stochastic/Stochastic_Stochastic.template"

namespace Sink = Source;

#endif
