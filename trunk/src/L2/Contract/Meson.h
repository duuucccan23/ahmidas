#pragma once

#include <cassert>

#include <L0/Base/Base.h>
#include <L0/Dirac/Gamma.h>
#include <L0/QCD/Tensor.h>
#include <L0/Core/Propagator.h>
#include <L0/Core/Correlator.h>


namespace Contract
{

  template< size_t IndexSrc, size_t IndexSnk >
  Core::Correlator light_meson_twopoint(Core::Propagator const *u, Core::Propagator const *d,
                                        Dirac::Gamma< IndexSrc > const &interpolSrc,
                                        Dirac::Gamma< IndexSnk > const &interpolSnk);

  template< size_t IndexSrc, size_t IndexSnk >
  Core::Correlator light_meson_twopoint(Core::Propagator const *u, Core::Propagator const *d,
                                        Dirac::Gamma< IndexSrc > const &interpolSrc,
                                        Dirac::Gamma< IndexSnk > const &interpolSnk,
                                        size_t const *momentum);
}

#include "Meson/Meson.inlines"
