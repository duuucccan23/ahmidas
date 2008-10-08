#ifndef GUARD_BASE_WEAVE_H
#define GUARD_BASE_WEAVE_H

#include <L0/Base/Base.h>

namespace Base
{
  template< size_t L, size_t T >
  class Weave
  {
    static Weave< L, T> *s_Weave;
    Weave< L, T>();

    public:
      static Weave< L, T > &instance();

      size_t localVolume() const;
      size_t dim(size_t idx) const;

      template< typename Element >
      void fieldShift(Base::SpaceTimeIndex idx, Base::Direction dir, Element *field, size_t const *offsets) const;
  };
}

#include "Weave/Weave.inlines"
#include "Weave/Weave.template"
#include "Weave/Weave_Weave.template"
#include "Weave/Weave_instance.template"

#endif
