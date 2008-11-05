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
      static Weave< L, T > &instance();

      size_t localVolume() const;
      size_t localSurface(Base::SpaceTimeIndex idx) const;
      size_t localSpatialVolume() const;
      size_t dim(Base::SpaceTimeIndex idx) const;
      size_t localSize(Base::SpaceTimeIndex idx) const;

      template< typename Element >
      void fieldShift(Base::SpaceTimeIndex idx, Base::Direction dir, Element *field, size_t const *offsets) const;

      size_t globalCoordToLocalIndex(size_t x, size_t y, size_t z) const;
      size_t globalCoordToLocalIndex(size_t x, size_t y, size_t z, size_t t) const;

      template< typename Element >
      void sumOverTimeSlices(Element *data) const;

    private:
      Weave< L, T>();
      bool isLocallyAvailable(size_t x, size_t y, size_t z) const;
      bool isLocallyAvailable(size_t x, size_t y, size_t z, size_t t) const;
      size_t fromGlobal(size_t x, Base::SpaceTimeIndex idx) const;
  };
}

#include "Weave/Weave.inlines"
#include "Weave/Weave.template"
#include "Weave/Weave_Weave.template"
#include "Weave/Weave_instance.template"

#endif
