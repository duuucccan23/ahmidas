#ifndef GUARD_SMEARING_FUZZ_H
#define GUARD_SMEARING_FUZZ_H

#include "Fields/GaugeField/GaugeField.h"
#include "Fields/GaugeComponent/GaugeComponent.h"
#include "SU3/Matrix/Matrix.h"

namespace Smearing
{
  class Fuzz
  {
    size_t  d_length;
    
    public:
      Fuzz(size_t length);
      void smear(Fields::GaugeField &field) const;
    
    private:
      void accumDirection(Fields::GaugeField &field, SpaceTimeIndex idx) const;
  };
}

#endif
