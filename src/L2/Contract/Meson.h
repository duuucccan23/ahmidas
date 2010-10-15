#pragma once

#include <vector>
#include <utility>
#include <L0/Dirac/Gamma.h>
#include <L0/Core/Propagator.h>
#include <L0/Core/Correlator.h>

namespace Contract
{
  template< size_t IndexSrc, size_t IndexSnk >
  Core::Correlator< Dirac::Matrix > light_meson_twopoint(Core::Propagator const *u, Core::Propagator const *d,
                                        Dirac::Gamma< IndexSrc > const &interpolSrc,
                                        Dirac::Gamma< IndexSnk > const &interpolSnk);

  template< size_t IndexSrc, size_t IndexSnk >
  Core::Correlator< Dirac::Matrix > light_meson_twopoint(Core::Propagator const *u, Core::Propagator const *d,
                                        Dirac::Gamma< IndexSrc > const &interpolSrc,
                                        Dirac::Gamma< IndexSnk > const &interpolSnk,
                                        size_t const *momentum);

  std::vector< Core::Correlator< Dirac::Matrix > > light_meson_twopoint(
    Core::Propagator const &psi1, Core::Propagator const &psi2,
    std::vector< std::pair< Base::Operator, Base::Operator > > const &operator_combinations);

  // this works using the one-end trick and gives all 16 gamma-combinations
  // return value: Array of 16 Core::Correlator arranged as follows:
  // gamma5, gamma0, gamma1, gamma2, gamma3, unity,
  // gamma5*gamma0, gamma5*gamma1, gamma5*gamma2, gamma5*gamma3,
  // sigma01, sigma02, sigma03, sigma12, sigma13, sigma23
  std::vector< Core::Correlator< Dirac::Matrix > > light_meson_twopoint_stochastic(
    Core::StochasticPropagator< 4 > const &psi1, Core::StochasticPropagator< 4 > const &psi2);

  std::vector< Core::Correlator< Dirac::Matrix > > light_meson_twopoint_stochastic(
    Core::StochasticPropagator< 4 > const &psi1, Core::StochasticPropagator< 4 > const &psi2,
    std::vector< std::pair< Base::Operator, Base::Operator > > const &operator_combinations);
}

#include "Meson/Meson.inlines"
