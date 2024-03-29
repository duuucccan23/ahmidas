#include "Disconnected.ih"


//Compute the disconnected loop of the operator O_{mumu} and O_{mumu}\gamma5
// O_00 and g5 O_00 cross checked independtly using an independent code 
// still to cross check the spatial derivative
//this function is used to perfom contraction on a  Core::StochasticPropagator< 1 >  
// not that it is very slow due to the absence of any iterator, and because the multiplication Core::StochasticPropagator< 1 >  is performed filling a Propagtors with 12 component.

namespace Contract
{
	
	// compute xi^* x Gamma x psi
	std::vector< Core::Correlator< Dirac::Matrix > > compute_loop_twist2_operator(
			Core::Field < QCD::Gauge > &gauge_field,
			Core::StochasticPropagator< 1 > const &xi, Core::StochasticPropagator< 1 > const &psi)
	{

		Dirac::Gamma<5> gamma5;
		Dirac::Gamma<1> gamma1;
		Dirac::Gamma<2> gamma2;
		Dirac::Gamma<3> gamma3;
		Dirac::Gamma<4> gamma4;

		size_t const L(psi.L());
		size_t const T(psi.T());

		std::vector< Core::Correlator< Dirac::Matrix > > twopoints;

		assert(xi.T() == psi.T() && psi.L() == psi.L());

		Core::StochasticPropagator< 1 > xi_conj(xi);
		xi_conj.conjugate();
		xi_conj.isolate();

		Core::StochasticPropagator< 1 > Gamma_psi(L,T);

		Core::StochasticPropagator< 1 > xi_conj_shifted_T_UP(xi_conj);
		Core::StochasticPropagator< 1 > xi_conj_shifted_T_DOWN(xi_conj);
		Core::StochasticPropagator< 1 > xi_conj_shifted_X_UP(xi_conj);
		Core::StochasticPropagator< 1 > xi_conj_shifted_X_DOWN(xi_conj);
		Core::StochasticPropagator< 1 > xi_conj_shifted_Y_UP(xi_conj);
		Core::StochasticPropagator< 1 > xi_conj_shifted_Y_DOWN(xi_conj);
		Core::StochasticPropagator< 1 > xi_conj_shifted_Z_UP(xi_conj);
		Core::StochasticPropagator< 1 > xi_conj_shifted_Z_DOWN(xi_conj);
		
		xi_conj_shifted_T_UP.isolate();
		xi_conj_shifted_T_DOWN.isolate();
		xi_conj_shifted_X_UP.isolate();
		xi_conj_shifted_X_DOWN.isolate();
		xi_conj_shifted_Y_UP.isolate();
		xi_conj_shifted_Y_DOWN.isolate();
		xi_conj_shifted_Z_UP.isolate();
		xi_conj_shifted_Z_DOWN.isolate();

		xi_conj_shifted_T_UP.shift(Base::idx_T, Base::dir_UP);
		xi_conj_shifted_T_DOWN.shift(Base::idx_T, Base::dir_DOWN);
		xi_conj_shifted_X_UP.shift(Base::idx_X, Base::dir_UP);
		xi_conj_shifted_X_DOWN.shift(Base::idx_X, Base::dir_DOWN);
		xi_conj_shifted_Y_UP.shift(Base::idx_Y, Base::dir_UP);
		xi_conj_shifted_Y_DOWN.shift(Base::idx_Y, Base::dir_DOWN);
		xi_conj_shifted_Z_UP.shift(Base::idx_Z, Base::dir_UP);
		xi_conj_shifted_Z_DOWN.shift(Base::idx_Z, Base::dir_DOWN);



		// now start the very ugly part ... O_mu mu and gamma_5 O_mu mu 

		//O44  

		{ // part1 tr g0 U_0(x) psi(x+0) xi^star(x)
			Core::Propagator tmp(psi);
			tmp.isolate();
			tmp.shift(Base::idx_T, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_T);
			tmp.rightMultiply(gamma4);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);
		}

		{ // part2 tr g0 U^dag_0(x-0) psi(x-0) xi^star(x)

			Core::Propagator tmp(psi);
			tmp.isolate();
			tmp.rightMultiplyDagger(gauge_field,Base::idx_T);  
			tmp.shift(Base::idx_T, Base::dir_UP);
			tmp.rightMultiply(gamma4);

			Gamma_psi= tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);
		}

		{// part3 tr g0 U^dag_0(x) psi(x) xi^star(x+0)
			Core::Propagator tmp(psi);
			
			tmp.isolate();
			tmp.rightMultiplyDagger(gauge_field,Base::idx_T); 
			tmp.rightMultiply(gamma4);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj_shifted_T_DOWN)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);

		}



		{// part4 tr g0 U_0(x-0) psi(x) xi^star(x-0)
			Core::Propagator tmp(psi);

			tmp.isolate();

			tmp.shift(Base::idx_T, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_T);
			tmp.rightMultiply(gamma4);
			tmp.shift(Base::idx_T, Base::dir_UP);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj_shifted_T_UP)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);


		}

		//end O44

		//O11
		{
			Core::Propagator tmp(psi);

			tmp.isolate();
			tmp.shift(Base::idx_X, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_X);
			tmp.rightMultiply(gamma1);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);
		}


		{
			Core::Propagator tmp(psi);

			tmp.isolate();
			tmp.rightMultiplyDagger(gauge_field,Base::idx_X);  
			tmp.shift(Base::idx_X, Base::dir_UP);
			tmp.rightMultiply(gamma1);

			Gamma_psi= tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);
		}

		{
			Core::Propagator tmp(psi);
			tmp.isolate();

			tmp.rightMultiplyDagger(gauge_field,Base::idx_X); 
			tmp.rightMultiply(gamma1);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj_shifted_X_DOWN)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);

		}

		{
			Core::Propagator tmp(psi);
			tmp.isolate();

			tmp.shift(Base::idx_X, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_X);
			tmp.rightMultiply(gamma1);
			tmp.shift(Base::idx_X, Base::dir_UP);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj_shifted_X_UP)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);


		}
		//end O11
		//O22
		{
			Core::Propagator tmp(psi);

			tmp.isolate();
			tmp.shift(Base::idx_Y, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_Y);
			tmp.rightMultiply(gamma2);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);
		}


		{
			Core::Propagator tmp(psi);

			tmp.isolate();
			tmp.rightMultiplyDagger(gauge_field,Base::idx_Y);  
			tmp.shift(Base::idx_Y, Base::dir_UP);
			tmp.rightMultiply(gamma2);

			Gamma_psi= tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);
		}

		{
			Core::Propagator tmp(psi);
			tmp.isolate();

			tmp.rightMultiplyDagger(gauge_field,Base::idx_Y); 
			tmp.rightMultiply(gamma2);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj_shifted_Y_DOWN)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);

		}

		{
			Core::Propagator tmp(psi);

			tmp.isolate();

			tmp.shift(Base::idx_Y, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_Y);
			tmp.rightMultiply(gamma2);
			tmp.shift(Base::idx_Y, Base::dir_UP);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj_shifted_Y_UP)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);


		}
		//end O22

		//O33 
		{
			Core::Propagator tmp(psi);

			tmp.isolate();
			tmp.shift(Base::idx_Z, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_Z);
			tmp.rightMultiply(gamma3);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);
		}


		{
			Core::Propagator tmp(psi);

			tmp.isolate();
			tmp.rightMultiplyDagger(gauge_field,Base::idx_Z);  
			tmp.shift(Base::idx_Z, Base::dir_UP);
			tmp.rightMultiply(gamma3);

			Gamma_psi= tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);
		}


		{
			Core::Propagator tmp(psi);
			tmp.isolate();

			tmp.rightMultiplyDagger(gauge_field,Base::idx_Z); 
			tmp.rightMultiply(gamma3);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj_shifted_Z_DOWN)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);

		}

		{
			Core::Propagator tmp(psi);
			tmp.isolate();

			tmp.shift(Base::idx_Z, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_Z);
			tmp.rightMultiply(gamma3);
			tmp.shift(Base::idx_Z, Base::dir_UP);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj_shifted_Z_UP)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);


		}
		//end O33



		//O44 polarized
		{
			Core::Propagator tmp(psi);
			tmp.isolate();
			tmp.shift(Base::idx_T, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_T);
			tmp.rightMultiply(gamma5);
			tmp.rightMultiply(gamma4);
			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);
		}



		{
			Core::Propagator tmp(psi);

			tmp.isolate();
			tmp.rightMultiply(gamma5); 
			tmp.rightMultiplyDagger(gauge_field,Base::idx_T);  
			tmp.shift(Base::idx_T, Base::dir_UP);
			tmp.rightMultiply(gamma4);

			Gamma_psi= tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);
		}

		{
			Core::Propagator tmp(psi);
			tmp.isolate();
		
			tmp.rightMultiply(gamma5);
			tmp.rightMultiplyDagger(gauge_field,Base::idx_T); 
			tmp.rightMultiply(gamma4);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj_shifted_T_DOWN)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);

		}

		{
			Core::Propagator tmp(psi);
			tmp.isolate();
			
			tmp.rightMultiply(gamma5);
			tmp.shift(Base::idx_T, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_T);
			tmp.rightMultiply(gamma4);
			tmp.shift(Base::idx_T, Base::dir_UP);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj_shifted_T_UP)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);


		}
		//end O44 polarized



		//O11 polarized
		{
			Core::Propagator tmp(psi);

			tmp.isolate();
			tmp.rightMultiply(gamma5);
			tmp.shift(Base::idx_X, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_X);
			tmp.rightMultiply(gamma1);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);
		}


		{
			Core::Propagator tmp(psi);

			tmp.isolate();
			tmp.rightMultiply(gamma5);
			tmp.rightMultiplyDagger(gauge_field,Base::idx_X);  
			tmp.shift(Base::idx_X, Base::dir_UP);
			tmp.rightMultiply(gamma1);

			Gamma_psi= tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);
		}

		{
			Core::Propagator tmp(psi);
			tmp.isolate();
			
			tmp.rightMultiply(gamma5);
			tmp.rightMultiplyDagger(gauge_field,Base::idx_X); 
			tmp.rightMultiply(gamma1);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj_shifted_X_DOWN)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);

		}

		{
			Core::Propagator tmp(psi);
			
			tmp.isolate();
		
			tmp.rightMultiply(gamma5);
			tmp.shift(Base::idx_X, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_X);
			tmp.rightMultiply(gamma1);
			tmp.shift(Base::idx_X, Base::dir_UP);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj_shifted_X_UP)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);


		}
		//end O11 polarized


		//O22 polarized
		{
			Core::Propagator tmp(psi);

			tmp.isolate();
			tmp.rightMultiply(gamma5);
			tmp.shift(Base::idx_Y, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_Y);
			tmp.rightMultiply(gamma2);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);
		}


		{
			Core::Propagator tmp(psi);

			tmp.isolate();
			tmp.rightMultiply(gamma5);
			tmp.rightMultiplyDagger(gauge_field,Base::idx_Y);  
			tmp.shift(Base::idx_Y, Base::dir_UP);
			tmp.rightMultiply(gamma2);

			Gamma_psi= tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);
		}

		{
			Core::Propagator tmp(psi);
			tmp.isolate();
			
			tmp.rightMultiply(gamma5);
			tmp.rightMultiplyDagger(gauge_field,Base::idx_Y); 
			tmp.rightMultiply(gamma2);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj_shifted_Y_DOWN)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);

		}

		{
			Core::Propagator tmp(psi);
			tmp.isolate();
			
			tmp.rightMultiply(gamma5);
			tmp.shift(Base::idx_Y, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_Y);
			tmp.rightMultiply(gamma2);
			tmp.shift(Base::idx_Y, Base::dir_UP);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj_shifted_Y_UP)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);


		}
		//end O22 polarized


		//O33 polarized
		{
			Core::Propagator tmp(psi);

			tmp.isolate();
			tmp.rightMultiply(gamma5);
			tmp.shift(Base::idx_Z, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_Z);
			tmp.rightMultiply(gamma3);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);
		}


		{
			Core::Propagator tmp(psi);

			tmp.isolate();
			tmp.rightMultiply(gamma5);
			tmp.rightMultiplyDagger(gauge_field,Base::idx_Z);  
			tmp.shift(Base::idx_Z, Base::dir_UP);
			tmp.rightMultiply(gamma3);

			Gamma_psi= tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);
		}

		{
			Core::Propagator tmp(psi);
			tmp.isolate();
		
			tmp.rightMultiply(gamma5);
			tmp.rightMultiplyDagger(gauge_field,Base::idx_Z); 
			tmp.rightMultiply(gamma3);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj_shifted_Z_DOWN)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);

		}

		{
			Core::Propagator tmp(psi);
			
			tmp.isolate();
			
			tmp.rightMultiply(gamma5);
			tmp.shift(Base::idx_Z, Base::dir_DOWN);
			tmp.rightMultiply(gauge_field,Base::idx_Z);
			tmp.rightMultiply(gamma3);
			tmp.shift(Base::idx_Z, Base::dir_UP);

			Gamma_psi = tmp;

			Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj_shifted_Z_UP)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));

			twopoint.sumOverSpatialVolume(); 
			twopoints.push_back(twopoint);


		}
		//end O33 polarized






		//contract

		return twopoints;
	}

}
