#infdef GUARD_SOURCE_POINT_H
#define GUARD_SOURCE_POINT_H

#include <L0/Base/Base.h>

namespace Source
{
  template< size_t L, size_t T >
  class Point
  {
    size_t d_coord[3];

    public:
      Point(size_t const x, size_t const y, size_t const z);
      coord(Base::SpaceTimeIndex const idx);
  };
}

#include "Point/Point.inlines"
#include "Point/Point_coord.template"

#endif
