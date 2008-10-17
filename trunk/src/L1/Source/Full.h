#ifndef GUARD_SOURCE_FULL_H
#define GUARD_SOURCE_FULL_H

#include <L0/Core/Field.h>
#include <L0/QCD/Spinor.h>

namespace Source
{
  template< size_t L, size_t T >
  class Full
  {
    Field< QCD::Spinor, L, T > d_source;

    public:
      Full(Field< QCD::Spinor, L, T > const &source);
  };
}

namespace Sink = Source;

#endif
