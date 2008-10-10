#ifndef GUARD_SMEAR_APE_H
#define GUARD_SMEAR_APE_H

#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/SU3/Matrix.h>
#include <L1/Path.h>
#include <L1/Tool.h>

// Performs an APE smearing step.
namespace Smear
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
