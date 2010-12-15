#include "Disconnected.ih"

namespace Contract
{

	// compute xi^* x Gamma x psi
	std::vector< Core::Correlator< Dirac::Matrix > > compute_loop_twist2_operator(
			Core::Field < QCD::Gauge > &gauge_field,
			Core::StochasticPropagator< 1 > const &xi, Core::StochasticPropagator< 1 > const &psi)
	{

		size_t const L(psi.L());
		size_t const T(psi.T());

		assert(xi.T() == psi.T() && psi.L() == psi.L());
		std::vector< Core::Correlator< Dirac::Matrix > > twopoints;

		Dirac::Gamma<5> gamma5;
		Dirac::Gamma<1> gamma1;
		Dirac::Gamma<2> gamma2;
		Dirac::Gamma<3> gamma3;
		Dirac::Gamma<4> gamma4;



		Core::StochasticPropagator< 1 > xi_conj(xi);
		//conjugation
		xi_conj.conjugate();


		Core::StochasticPropagator< 1 > Gamma_psi(L,T);
		Core::StochasticPropagator<1>   prop(xi_conj);
		// now start the very ugly part ... O_mu mu and gamma_5 O_mu mu 

		//O44 

		{
			Core::Propagator tmp(psi);
			tmp.isolate();
			tmp.shift(Base::idx_T, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_T);
			tmp.rightMultiply(gamma4);
			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			//sum over space
			twopoint.sumOverSpatialVolume(); 
			//accumulate 
			twopoints.push_back(twopoint);
		}


		{
			Core::Propagator tmp(psi);
			tmp.isolate();
			tmp.rightMultiplyDagger(gauge_field,Base::idx_T);  
			tmp.shift(Base::idx_T, Base::dir_UP);
			tmp.rightMultiply(gamma4);
			Gamma_psi= tmp;
			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));
			//sum over space
			twopoint.sumOverSpatialVolume(); 
			//accumulate 
			twopoints.push_back(twopoint);
		}

		{
			Core::Propagator tmp(Gamma_psi);

			prop.isolate();
			tmp.isolate();
			tmp.rightMultiplyDagger(gauge_field,Base::idx_T); 
		//	tmp.shift(Base::idx_T, Base::dir_UP);
			tmp.rightMultiply(gamma4);
			prop.shift(Base::idx_T, Base::dir_UP);
			Gamma_psi =tmp;
			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(prop)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));
			//sum over space
			twopoint.sumOverSpatialVolume(); 
			//accumulate 
			twopoints.push_back(twopoint);

		}

		{
			Core::Propagator tmp(prop);
			tmp.isolate();
		//	tmp.shift(Base::idx_T, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_T);
			tmp.rightMultiply(gamma4);
			prop.shift(Base::idx_T, Base::dir_DOWN);
			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(prop)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));
			//sum over space
			twopoint.sumOverSpatialVolume(); 
			//accumulate 
			twopoints.push_back(twopoint);



		}





		//contract

		return twopoints;
	}

}
