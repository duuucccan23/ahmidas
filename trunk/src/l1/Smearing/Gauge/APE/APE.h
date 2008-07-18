#ifndef GUARD_SMEARING_APE_H
#define GUARD_SMEARING_APE_H

#include <l0/Core/Buffer/Buffer.h>
#include <l0/Core/Field/Field.h>
#include <l0/Core/Component/Component.h>
#include <l0/QCD/Gauge/Gauge.h>
#include <l0/SU3/Matrix/Matrix.h>

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

#include "APE.inlines"
#include "APE_a.template"
#include "APE_b.template"

#endif
