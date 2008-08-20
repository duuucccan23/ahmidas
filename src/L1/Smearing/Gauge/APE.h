#ifndef GUARD_SMEARING_APE_H
#define GUARD_SMEARING_APE_H

#include <L0/Core/Buffer.h>
#include <L0/Core/Field.h>
#include <L0/Core/Component.h>
#include <L0/QCD/Gauge.h>
#include <L0/SU3/Matrix.h>

// Performs an APE smearing step.
namespace Smearing
{
  class APE
  {
    double d_alpha;

    public:
      APE(double alpha);

      template< size_t L, size_t T >
      void smear(Core::Field< QCD::Gauge, L, T > &field) const;

      template< size_t L, size_t T >
      void smear(Core::Field< QCD::Gauge, L, T > &field, size_t iterations) const;
  };
}

#include "APE/APE.inlines"
#include "APE/APE_smear_a.template"
#include "APE/APE_smear_b.template"

#endif
