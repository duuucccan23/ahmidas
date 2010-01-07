#pragma once

#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>

// Performs an APE smearing step.
namespace Smear
{
  class APE
  {
    double d_alpha;

    public:
      APE(double alpha);

      void smear(Core::Field< QCD::Gauge > &field) const;
      void smear(Core::Field< QCD::Gauge > &field, size_t iterations) const;
  };
}

#include "APE/APE.inlines"
