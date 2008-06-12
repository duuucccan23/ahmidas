#ifndef GUARD_SMEARING_APE_H
#define GUARD_SMEARING_APE_H

#include "Fields/GaugeField/GaugeField.h"
#include "Fields/GaugeComponent/GaugeComponent.h"
#include "SU3/Matrix/Matrix.h"

// Performs an APE smearing step.
namespace Smearing
{
  class APE
  {
    double d_alpha;
    
    public:
      APE(double alpha);
      void smear(Fields::GaugeField &field) const;
      void smear(Fields::GaugeField &field, size_t iterations) const;
  };
}

#include "APE.inlines"

#endif
