#ifndef GUARD_PATH_H
#define GUARD_PATH_H

#include <L0/Base/Base.h>
#include <L0/Core/Component.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/SU3/Matrix.h>

namespace Path
{
  template< size_t L, size_t T >
  Core::Field< SU3::Matrix, L, T > staple(Core::Field< QCD::Gauge, L, T > &field, Base::SpaceTimeIndex via, Base::Direction dirVia, Base::SpaceTimeIndex to, Base::Direction dirTo);

  template< size_t L, size_t T >
  Core::Field< SU3::Matrix, L, T > square(Core::Field< QCD::Gauge, L, T > &field, Base::SpaceTimeIndex via, Base::Direction dirVia, Base::SpaceTimeIndex to, Base::Direction dirTo);
}

#include "Path/Path_staple.template"
#include "Path/Path_square.template"

#endif

