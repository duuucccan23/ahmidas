#ifndef GUARD_TOOL_H
#define GUARD_TOOL_H

#include <L0/Base/Weave.h>
#include <L0/Core/Field.h>
#include <L0/Core/Propagator.h>
#include <L0/QCD/Gauge.h>
#include <L0/SU3/Matrix.h>

namespace Tool
{
  template< size_t L, size_t T >
  void reunitarize(Core::Field< QCD::Gauge, L, T > *field);

  template< typename Element, size_t L, size_t T >
  void randomize(Core::Field< Element, L, T > *field);

  template< size_t L, size_t T >
  void randomize(Core::Propagator< L, T > *propagator);
  
  template< typename Element, size_t L, size_t T >
  void setToZero(Core::Field< Element, L, T > *field);

}

#include "Tool/Tool.inlines"

#include "Tool/Tool_reunitarize.template"
#include "Tool/Tool_randomize.template"
#include "Tool/Tool_setToZero.template"

#endif
