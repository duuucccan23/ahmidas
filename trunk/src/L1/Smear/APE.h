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
      void smear(Core::Field< QCD::Gauge > &field, size_t const iterations) const;

      // FIX THIS: those routine is just a workaround for better performance in scalar code
      void smear(Core::Field< QCD::Gauge > &field, size_t const iterations, size_t const timeslice) const;
  };
}

#include "APE/APE.inlines"
