#ifndef GUARD_SMEARING_APE_H
#define GUARD_SMEARING_APE_H

#include <Core/Buffer/Buffer.h>
#include <Core/Field/Field.h>
#include <Core/Component/Component.h>
#include <QCD/Gauge/Gauge.h>
#include <SU3/Matrix/Matrix.h>

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
