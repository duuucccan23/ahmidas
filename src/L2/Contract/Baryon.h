#pragma once

#include <vector>
#include <L0/QCD/Gauge.h>
#include <L0/Core/Field.h>
#include <L0/Core/Propagator.h>
#include <L0/Core/Correlator.h>

namespace Contract
{

  Core::Correlator proton_twopoint(Core::Propagator const &u1, Core::Propagator const &u2, Core::Propagator const &d,
                                   Base::BaryonPropagatorProjector const projector);


  void create_sequential_source_proton_d(Core:: Propagator * const seqSrc,
                                         Core::Propagator const &u1, Core::Propagator const &u2, size_t const t_snk);

  // the proton three point with a stochastic (estimate of an) all-to-all propagator at the proton sink.
  // this is a mess concerning performance since it scales vith 4-Volume^2, but provides a good cross-check
  // the gauge field is sometimes not needed, in this case NULL can be passed;
  std::vector< Core::Correlator > proton_threepoint_stochastic_naive(Core::Propagator const &u,
                                                          Core::Propagator const &d,
                                                          Core::StochasticPropagator<12> const &u_stoch_at_sink,
                                                          Core::StochasticPropagator<12> const &d_stoch_at_sink,
                                                          Core::StochasticSource<12> const &xi_at_sink,
                                                          Core::Field < QCD::Gauge > const * const gauge_field,
                                                          std::vector< Base::Operator > const &ops,
                                                          size_t const t_src, size_t const t_snk);

  std::vector< Core::Correlator > proton_threepoint_sequential(
    Core:: Propagator const &bw_prop_u, Core::Propagator const &fw_prop_u,
    Core:: Propagator const &bw_prop_d, Core::Propagator const &fw_prop_d,
    std::vector< Base::Operator > ops, Base::BaryonPropagatorProjector const my_projector);

  std::vector< Core::Correlator > proton_threepoint_stochastic(Core::Propagator const &u,
                                                Core::Propagator const &d,
                                                Core::StochasticPropagator <12> const &phi_u,
                                                Core::StochasticPropagator <12> const &phi_d,
                                                Core::StochasticSource <12> const &xi,
                                                size_t t_source, size_t t_sink,
                                                /* this allows for more than one operator */
                                                std::vector< Base::Operator > const &my_operators,
                                                Base::BaryonPropagatorProjector const my_projector);
}
