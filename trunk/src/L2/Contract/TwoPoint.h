#pragma once

#include <L0/Core/Correlator.h>
#include <L0/Core/TMatrix.h>

namespace Contract
{
  template< size_t L, size_t T >
  Core::Correlator< L, T > twoPoint(Core::hcTMatrix< L, T > const &left, Core::TMatrix< L, T > const &right);
}

#include "TwoPoint/twoPoint_a.template"
