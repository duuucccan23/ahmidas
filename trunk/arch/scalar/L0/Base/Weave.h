#pragma once

#include<complex>

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

      size_t globalCoordToLocalIndex(size_t const x, size_t const y, size_t const z) const;
      size_t globalCoordToLocalIndex(size_t const x, size_t const y, size_t const z, size_t const t) const;

      // This funtion exists only for compatibility with parallel implementation
      // here, it simply does not do anything but copy the data;
      void sumOverTimeSlices(std::complex< double > const *data_send, 
                             std::complex< double > *data_recv, size_t const count) const;

      bool isLocallyAvailable(size_t const x, size_t const y, size_t const z) const;
      bool isLocallyAvailable(size_t const x, size_t const y, size_t const z, size_t const t) const;

    private:
      size_t fromGlobal(size_t const idx, Base::SpaceTimeIndex const mu) const;
  };
}

#include "Weave/Weave.inlines"
