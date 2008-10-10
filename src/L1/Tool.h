#ifndef GUARD_TOOL_H
#define GUARD_TOOL_H

#include <L0/Base/Weave.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Gauge.h>
#include <L0/SU3/Matrix.h>

namespace Tool
{
  template< size_t L, size_t T >
  void reunitarize(Core::Field< QCD::Gauge, L, T > *field);

  template< size_t L, size_t T >
  void randomGauge(Core::Field< QCD::Gauge, L, T > *field);
}

#include "Tool/reunitarize.template"
#include "Tool/randomGauge.template"

#endif
