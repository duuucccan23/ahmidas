#pragma once

#include <L0/Core/Field.h>
#include <L0/Core/Component.h>
#include <L0/QCD/Gauge.h>
#include <L0/QCD/Spinor.h>
#include <L0/QCD/Tensor.h>
#include <L1/Transport.h>

namespace Smear
{
  class Jacobi
  {
    double d_kappa;
    double d_weight;

    public:
      Jacobi(double kappa);

      void smear(Core::Field< QCD::Spinor > *spinorField, Core::Field< QCD::Gauge > &gaugeField) const;

      void smear(Core::Field< QCD::Spinor > *spinorField, Core::Field< QCD::Gauge > &gaugeField, size_t const iterations) const;

      void smear(Core::Field< QCD::Tensor > *tensorField, Core::Field< QCD::Gauge > &gaugeField) const;

      void smear(Core::Field< QCD::Tensor > *tensorField, Core::Field< QCD::Gauge > &gaugeField, size_t const iterations) const;

      // FIX THIS: those routine is just a workaround for better performance in scalar code
      void smear(Core::Field< QCD::Tensor > *tensorField, Core::Field< QCD::Gauge > &gaugeField, size_t const iterations, size_t const timeslice) const;
      void smear(Core::Field< QCD::Spinor > *spinorField, Core::Field< QCD::Gauge > &gaugeField, size_t const iterations, size_t const timeslice) const;
  
//       void smear(Source::Point< L, T > *source, Core::Field< QCD::Gauge > &gaugeField, Base::ColourIndex, Base::DiracIndex) const;
// 
//       void smear(Source::Point< L, T > *source, Core::Field< QCD::Gauge > &gaugeField,
//                  Base::ColourIndex, Base::DiracIndex, size_t iterations) const;

  };
}

#include "Jacobi/Jacobi.inlines"
