#ifndef GUARD_CONTRACT_WICK_H
#define GUARD_CONTRACT_WICK_H

#include <L0/Core/Propagator.h>
#include <L0/Core/TMatrix.h>
#include <L0/QCD/Spinor.h>

namespace Contract
{
  template< size_t L, size_t T >
  Core::Correlator &wick(Core::Source::Point< L, T > source, Core::Field< QCD::Spinor, L, T > field);

  
  template< size_t L, size_t T >
  Core::TMatrix< T > &wick(Core::Source::Point< L, T >, Core::Propagator< L, T >);
}

#endif
