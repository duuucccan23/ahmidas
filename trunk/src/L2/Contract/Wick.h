#ifndef GUARD_CONTRACT_WICK_H
#define GUARD_CONTRACT_WICK_H

#include <L0/Core/Correlator.h>
#include <L0/Core/Propagator.h>
#include <L0/Core/TMatrix.h>
#include <L0/QCD/Spinor.h>
#include <L1/Source/Point.h>

namespace Contract
{
  template< size_t L, size_t T >
  Core::Correlator< T > &wick(Core::Source::Point< L, T > const &source, Base::DiracIndex const &dir,Base::ColourIndex const &col,
                              Core::Field< QCD::Spinor, L, T > const &field);

  template< size_t L, size_t T >
  Core::TMatrix< T > &wick(Core::Source::Point< L, T > const &source, Core::Propagator< L, T > const &propagator);
}

#include "Wick/wick_a.template"

#endif
