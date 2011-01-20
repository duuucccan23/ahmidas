#pragma once

#include <vector>
#include <utility>
#include <L0/Dirac/Gamma.h>
#include <L0/Core/Propagator.h>
#include <L0/Core/Correlator.h>
#include <L0/Core/Component.h>

namespace Contract
{
	std::vector< Core::Correlator< Dirac::Matrix > > compute_loop(
			Core::StochasticPropagator< 1 > const &xi, Core::StochasticPropagator< 1 > const &psi,
			std::vector< Base::Operator > ops);


	std::vector< Core::Correlator< Dirac::Matrix > > compute_loop(
			Core::StochasticPropagator< 1 > const &xi, Core::StochasticPropagator< 1 > const &psi,
			std::vector< Base::HermitianBilinearOperator > ops);

	std::vector< std::complex<double>  > compute_loop(
			Core::Propagator const &xi, Core::Propagator const &psi,
			std::vector< Base::HermitianBilinearOperator > ops);

	std::vector< std::complex<double>  > compute_loop_new(
			Core::Propagator const &xi, Core::Propagator const &psi,
			std::vector< Base::HermitianBilinearOperator > ops);

	std::vector< Core::Correlator< Dirac::Matrix > > compute_loop_twist2_operator(
			Core::Field < QCD::Gauge > &gauge_field,	
			Core::StochasticPropagator< 1 > const &xi, Core::StochasticPropagator< 1 > const &psi );

	std::vector< std::complex<double>  > compute_loop_twist2_operator(
			Core::Field < QCD::Gauge > &gauge_field,
			Core::Propagator const &xi, Core::Propagator const &psi);

}
