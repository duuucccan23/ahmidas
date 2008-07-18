#ifndef GUARD_SMEARING_FUZZ_H
#define GUARD_SMEARING_FUZZ_H

#include <algorithm>

#include <Core/Buffer/Buffer.h>
#include <Core/Field/Field.h>
#include <Core/Component/Component.h>
#include <QCD/Gauge/Gauge.h>
#include <SU3/Matrix/Matrix.h>

namespace Smearing
{
  class Fuzz
  {
    size_t  d_length;

    public:
      Fuzz(size_t length);

      template< size_t L, size_t T >
      void smear(Core::Field< QCD::Gauge, L, T > &field) const;

    private:
      template< size_t L, size_t T >      
      void accumDirection(Core::Field< QCD::Gauge, L, T > &field, Core::SpaceTimeIndex) const;
  };
}

#include "Fuzz.inlines"
#include "Fuzz_a.template"
#include "Fuzz_b.template"

#endif
