#ifndef GUARD_SOURCE_FULL_H
#define GUARD_SOURCE_FULL_H

#include <L0/Core/Field.h>
#include <L0/QCD/Spinor.h>

namespace Source
{
  template< size_t L, size_t T >
  class Full
  {
    Core::Field< QCD::Spinor, L, T > d_source;

    public:
      Full();
      Full(Core::Field< QCD::Spinor, L, T > const &source);

      Core::Field< QCD::Spinor, L, T > &source();
  };
}

#include "Full/Full.inlines"

namespace Sink = Source;

#endif
