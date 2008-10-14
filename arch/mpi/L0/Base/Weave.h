#ifndef GUARD_BASE_WEAVE_H
#define GUARD_BASE_WEAVE_H

#include <L0/Base/Base.h>
#include <L0/Base/Grid.h>
// This particular version of Weave is tailored for an MPI setup.

namespace Base
{
  template< size_t L, size_t T >
  class Weave
  {
    static Weave< L, T> *s_Weave;

    size_t d_surfaces[4]; //Surface size in direction
    size_t d_localVolume; //Local volume
    size_t d_localSize[4]; //Local size in direction

    public:
      Grid< L, T > d_grid;

      static Weave< L, T > &instance();
      size_t localVolume() const;
      size_t localSurface(Base::SpaceTimeIndex idx) const;
      size_t dim(Base::SpaceTimeIndex idx) const;
      size_t localSize(Base::SpaceTimeIndex idx) const;

      template< typename Element >
      void fieldShift(Base::SpaceTimeIndex idx, Base::Direction dir, Element *field, size_t const *offsets) const;

    private:
      Weave< L, T>();
  };
}

#include "Weave/Weave.template"
#include "Weave/Weave.inlines"
#include "Weave/Weave_instance.template"
#include "Weave/Weave_Weave.template"
#include "Weave/Weave_fieldShift.template"

#endif
