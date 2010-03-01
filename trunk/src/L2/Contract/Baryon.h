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

  Core::Correlator proton_twopoint(Core::Propagator const &u, Core::Propagator const &d,
                                   Base::BaryonPropagatorProjector const projector);
}
