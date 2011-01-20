#include "Disconnected.ih"

namespace Contract
{

  // compute xi^* x Gamma x psi
  std::vector< Core::Correlator< Dirac::Matrix > > compute_loop(
    Core::StochasticPropagator< 1 > const &xi, Core::StochasticPropagator< 1 > const &psi,
    std::vector< Base::HermitianBilinearOperator > ops)
  {

	

    assert(xi.T() == psi.T() && psi.L() == psi.L());

    std::vector< Core::Correlator< Dirac::Matrix > > twopoints;

	Core::StochasticPropagator< 1 > xi_conj(xi);
	
    //conjugation
	xi_conj.conjugate();

	// loop over all operator combinations

	for (size_t iOp=0; iOp<ops.size(); iOp++)
	{

		Core::StochasticPropagator< 1 > Gamma_psi(psi);

		//apply operator X
		Gamma_psi.rightMultiplyOperator(ops[iOp]);

		//contract
		Core::Correlator< Dirac::Matrix >twopoint( (dynamic_cast< Core::StochasticPropagator<1> & >(xi_conj)) * (dynamic_cast< Core::StochasticPropagator<1> & >(Gamma_psi)));
		//sum over space
		twopoint.sumOverSpatialVolume(); 


		twopoint.deleteField();	
		//accumulate 
		twopoints.push_back(twopoint);

		
	}

	return twopoints;
  }

}
