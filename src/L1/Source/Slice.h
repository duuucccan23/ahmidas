#pragma once

#include <L0/Core/Field.h>
#include <L0/QCD/Spinor.h>

namespace Source
{
  template< size_t L, size_t T >
  class Slice
  {
    Core::Field< QCD::Spinor, L, 1 >  d_source;

    public:
      Slice(Stochastic const &source, Base::ColourIndex, Base::DiracIndex);
  };
}

namespace Sink = Source;
