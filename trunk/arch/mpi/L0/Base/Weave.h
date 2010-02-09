#pragma once

#include<complex>

#include <L0/Base/Base.h>
#include <L0/Base/Grid.h>


// This particular version of Weave is tailored for an MPI setup.
namespace Base
{
  class Weave
  {
    size_t d_surfaces[4];
    size_t d_localSize[4];
    size_t d_L;
    size_t d_T;
    size_t d_localVolume;
    size_t d_globalVolume;

    public:
      Grid d_grid;

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

//       template< typename Element >
//       void sumOverTimeSlices(Element const *data_send, Element *data_recv);

      void sumOverTimeSlices(std::complex< double > const *data_send, std::complex< double > *data_recv)  const;

      // SD: Compiler complained about the const declaration, so I left it out
      // this could be cured by a const_cast of d_grid within the function body
      bool isLocallyAvailable(size_t x, size_t y, size_t z) const;
      bool isLocallyAvailable(size_t x, size_t y, size_t z, size_t t) const;

    private:
      size_t fromGlobal(size_t const x, Base::SpaceTimeIndex const idx) const;
  };
}

#include "Weave/Weave.inlines"
#include "Weave/Weave_fieldShift.template"
