#ifndef GUARD_SOURCE_BOX_H
#define GUARD_SOURCE_BOX_H

namespace Source
{
  template< size_t L, size_t T >
  class Box
  {
    size_t d_lower[3];
    size_t d_upper[3];
  };
}

namespace Sink = Source;

#endif
