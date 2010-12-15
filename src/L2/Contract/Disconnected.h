#pragma once

#include <vector>
#include <utility>
#include <L0/Dirac/Gamma.h>
#include <L0/Core/Propagator.h>
#include <L0/Core/Correlator.h>

namespace Contract
{
  std::vector< Core::Correlator< Dirac::Matrix > > compute_loop(
    Core::StochasticPropagator< 1 > const &xi, Core::StochasticPropagator< 1 > const &psi,
     std::vector< Base::Operator > ops);
  
  
  std::vector< Core::Correlator< Dirac::Matrix > > compute_loop(
    Core::StochasticPropagator< 1 > const &xi, Core::StochasticPropagator< 1 > const &psi,
     std::vector< Base::HermitianBilinearOperator > ops);

  
  std::vector< Core::Correlator< Dirac::Matrix > > compute_loop_twist2_operator(
		  Core::Field < QCD::Gauge > &gauge_field,	
		  Core::StochasticPropagator< 1 > const &xi, Core::StochasticPropagator< 1 > const &psi );


}
