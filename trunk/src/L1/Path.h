#pragma once

#include <L0/Base/Base.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/SU3/Matrix.h>

namespace Path
{
  void step(Core::Field< SU3::Matrix > *path, Core::Field< QCD::Gauge > &field, Base::SpaceTimeIndex idx, Base::Direction dir, size_t nsteps = 1);

  Core::Field< SU3::Matrix > step(Core::Field< QCD::Gauge > &field, Base::SpaceTimeIndex idx, Base::Direction dir, size_t nsteps = 1);

  Core::Field< SU3::Matrix > staple(Core::Field< QCD::Gauge > &field, Base::SpaceTimeIndex via, Base::Direction dirVia, Base::SpaceTimeIndex to, Base::Direction dirTo);

  Core::Field< SU3::Matrix > square(Core::Field< QCD::Gauge > &field, Base::SpaceTimeIndex via, Base::Direction dirVia, Base::SpaceTimeIndex to, Base::Direction dirTo);
}
