#pragma once

#include <L0/Base/Base.h>

namespace Base
{
  class Weave
  {

    size_t d_surfaces[4];
    size_t d_L;
    size_t d_T;
    size_t d_localVolume;
    size_t d_globalVolume;

    public:
      Weave(size_t const L, size_t const T);
      Weave(Weave const &other);
      Weave &operator=(Weave const &other);

      size_t L() const;
      size_t T() const;
      size_t localVolume() const;
      size_t localSurface(Base::SpaceTimeIndex idx) const;
      size_t localSpatialVolume() const;
      size_t dim(Base::SpaceTimeIndex idx) const;
      size_t localSize(Base::SpaceTimeIndex idx) const;
      size_t globalVolume() const;

      double sum(double result) const;

      template< typename Element >
      void fieldShift(Base::SpaceTimeIndex idx, Base::Direction dir, Element *field, size_t const *offsets) const;

      size_t globalCoordToLocalIndex(size_t x, size_t y, size_t z) const;
      size_t globalCoordToLocalIndex(size_t x, size_t y, size_t z, size_t t) const;

      template< typename Element >
      void sumOverTimeSlices(Element *data) const;

      bool isLocallyAvailable(size_t x, size_t y, size_t z) const;
      bool isLocallyAvailable(size_t x, size_t y, size_t z, size_t t) const;

    private:
      size_t fromGlobal(size_t x, Base::SpaceTimeIndex idx) const;
  };
}

#include "Weave/Weave.inlines"
