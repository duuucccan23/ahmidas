#pragma once

#include <L0/Base/Weave.h>
#include <L0/SU3/Vector.h>
#include <L0/Base/Random.h>
#include <cmath>

namespace Source
{
  template< size_t L, size_t T >
  class Stochastic
  {
    SU3::Vector *d_source;

    public:
      Stochastic();
      ~Stochastic();

      SU3::Vector const *source() const;
  };
}

#include "Stochastic/Stochastic.inlines"
#include "Stochastic/Stochastic_Stochastic.template"

namespace Sink = Source;
