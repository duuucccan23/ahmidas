#ifndef GUARD_BASE_WEAVE_H
#define GUARD_BASE_WEAVE_H

#include <L0/Base/Base.h>

namespace Base
{
  template< size_t L, size_t T >
  class Weave
  {
    static Weave< L, T> *s_Weave;
    size_t d_surfaces[4];

    public:
      bool const d_bigEndian;
      static Weave< L, T > &instance();

      size_t localVolume() const;
      size_t localSurface(Base::SpaceTimeIndex idx) const;
      size_t dim(Base::SpaceTimeIndex idx) const;
      size_t localSize(Base::SpaceTimeIndex idx) const;

      template< typename Element >
      void fieldShift(Base::SpaceTimeIndex idx, Base::Direction dir, Element *field, size_t const *offsets) const;

      size_t globalCoordToLocalIndex(size_t * const global) const;

      template< typename Element >
      void sumOverTimeSlices(Element *data) const;

    private:
      Weave< L, T>();
  };
}

#include "Weave/Weave.inlines"
#include "Weave/Weave.template"
#include "Weave/Weave_Weave.template"
#include "Weave/Weave_instance.template"

#endif
