#ifndef GUARD_SMEARING_JACOBI_H
#define GUARD_SMEARING_JACOBI_H

#include <L0/Core/Buffer.h>
#include <L0/Core/Field.h>
#include <L0/Core/Component.h>
#include <L0/QCD/Gauge.h>
#include <L0/QCD/Spinor.h>

namespace Smearing
{
  class Jacobi
  {
    double d_kappa;
    double d_weight;
    
    public:
      Jacobi(double kappa);
      
      template< size_t L, size_t T >
      void smear(Core::Field< QCD::Spinor, L, T > &spinorField, Core::Field< QCD::Gauge, L, T > &gaugeField) const;

      template< size_t L, size_t T >
      void smear(Core::Field< QCD::Spinor, L, T > &spinorField, Core::Field< QCD::Gauge, L, T > &gaugeField, size_t iterations) const;
  };
}

#include "Jacobi/Jacobi.inlines"
#include "Jacobi/Jacobi_smear.template"

#endif
