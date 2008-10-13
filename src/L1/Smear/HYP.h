#ifndef GUARD_SMEAR_HYP_H
#define GUARD_SMEAR_HYP_H

#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/SU3/Matrix.h>
#include <L1/Path.h>
#include <L1/Tool.h>

// Performs an HYP smearing step.
namespace Smear
{
  class HYP
  {
    double d_alpha;

    public:
      HYP(double alpha);

      template< size_t L, size_t T >
      void smear(Core::Field< QCD::Gauge, L, T > &field) const;

      template< size_t L, size_t T >
      void smear(Core::Field< QCD::Gauge, L, T > &field, size_t iterations) const;
  };
}

#include "HYP/HYP.inlines"
#include "HYP/HYP_smear_a.template"
#include "HYP/HYP_smear_b.template"

#endif
