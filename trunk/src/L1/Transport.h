#pragma once

#include <L0/Base/Base.h>
#include <L0/Core/Component.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/QCD/Spinor.h>
#include <L0/SU3/Matrix.h>

namespace Transport
{
  template< size_t L, size_t T >
  Core::Field< QCD::Spinor, L, T > step(Core::Field< QCD::Spinor, L, T > const &spinor,
      Core::Field< QCD::Gauge, L, T > const &gauge,
      Base::SpaceTimeIndex idx, Base::Direction dir);

  template< size_t L, size_t T >
  Core::Field< QCD::Spinor, L, T > range(Core::Field< QCD::Spinor, L, T > const &spinor,
      Core::Field< QCD::Gauge, L, T > const &gauge,
      Base::SpaceTimeIndex idx, Base::Direction dir, size_t steps);
}

#include "Transport_step.template"
#include "Transport_range.template"
