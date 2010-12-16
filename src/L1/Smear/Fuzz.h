#pragma once

#include <algorithm>

#include <L0/Core/Component.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/SU3/Matrix.h>

namespace Smear
{
  class Fuzz
  {
    size_t  d_length;

    public:
      Fuzz(size_t length);
      void smear(Core::Field< QCD::Gauge > &field) const;

    private:
      void accumDirection(Core::Field< QCD::Gauge > &field, Base::SpaceTimeIndex) const;
  };
}

#include "Fuzz/Fuzz.inlines"
//#include "Fuzz/Fuzz_accumDirection.template"
//#include "Fuzz/Fuzz_smear.template"
