#ifndef GUARD_SOURCE_STOCHASTIC_H
#define GUARD_SOURCE_STOCHASTIC_H

#include <L0/Core/Field.h>

namespace Source
{
  template< size_t L, size_t T >
  class Stochastic
  {
    Field< std::complex< double >, L, 1 >  d_source;

    public:
      Stochastic();
  };
}

#endif
