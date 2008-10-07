#ifndef GUARD_PATH_H
#define GUARD_PATH_H

#include <L0/Core/Component.h>
#include <L0/Core/Core.h>
#include <L0/Core/Buffer.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/SU3/Matrix.h>

namespace Path
{
  template< size_t L, size_t T >
  Core::Buffer< SU3::Matrix > staple(Core::Field< QCD::Gauge, L, T > &field, Core::SpaceTimeIndex towards, Core::Direction dirTo, Core::SpaceTimeIndex over, Core::Direction dirOver);

  template< size_t L, size_t T >
  Core::Buffer< std::complex< double > > plaquette(Core::Field< QCD::Gauge, L, T > &field, Core::SpaceTimeIndex towards, Core::Direction dirTo, Core::SpaceTimeIndex over, Core::Direction dirOver));
}

#include "Path/Path_staple.template"
#include "Path/Path_square.template"

#endif

