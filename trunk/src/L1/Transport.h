#pragma once

#include <L0/Base/Base.h>
#include <L0/Core/Component.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/QCD/Spinor.h>
#include <L0/SU3/Matrix.h>

namespace Transport
{

  Core::Field< QCD::Spinor > step(Core::Field< QCD::Spinor > const &spinor,
      Core::Field< QCD::Gauge > &gauge,
      Base::SpaceTimeIndex const idx, Base::Direction const dir);

  Core::Field< QCD::Spinor > range(Core::Field< QCD::Spinor > const &spinor,
      Core::Field< QCD::Gauge > &gauge,
      Base::SpaceTimeIndex const idx, Base::Direction const dir, size_t const steps);
}

#include "Transport/Transport_step.template"
#include "Transport/Transport_range.template"
