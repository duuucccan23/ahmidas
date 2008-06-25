#ifndef GUARD_SMEARING_JACOBI_H
#define GUARD_SMEARING_JACOBI_H

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
  }
}

#include "Jacobi.inlines"
#include "Jacobi_a.template"

#endif
