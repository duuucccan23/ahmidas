#pragma once

#include <cassert>
#include <string>
#include <vector>

#include <L0/Base/Base.h>
#include <L0/Dirac/Gamma.h>
#include <L0/QCD/Tensor.h>
#include <L0/Core/Propagator.h>
#include <L0/Core/Correlator.h>


namespace Contract
{

  Core::Correlator proton_twopoint(Core::Propagator const &u1, Core::Propagator const &u2, Core::Propagator const &d,
                                   Base::BaryonPropagatorProjector const projector);

  std::vector< Core::Correlator > proton_threepoint_stochastic(Core::Propagator const &u,
                                                Core::Propagator const &d,
                                                Core::StochasticPropagator <12> const &phi_u,
                                                Core::StochasticPropagator <12> const &phi_d,
                                                Core::StochasticSource <12> const &xi,
                                                size_t t_source, size_t t_sink,
                                                /* one eventually might skip this and iterate over all operators */
                                                Base::Operator my_operator,
                                                Base::BaryonPropagatorProjector const my_projector);
}
