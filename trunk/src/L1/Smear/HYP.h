#pragma once

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
    double d_alpha1;
    double d_alpha2;
    double d_alpha3;

    void smear(Core::Field< QCD::Gauge > &field) const;

    public:
      HYP(double alpha1, double alpha2, double alpha3);

      void smear(Core::Field< QCD::Gauge > &field, size_t iterations) const;
  };
}

#include "HYP/HYP.inlines"
