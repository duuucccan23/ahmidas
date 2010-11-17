#include "Meson.ih"

namespace Contract
{

  // this works using the one-end trick and gives all 16 gamma-combinations
  std::vector< Core::Correlator< Dirac::Matrix > > BK_threepoint_connected_stochastic(
    Core::StochasticPropagator< 4 > const &prop_dR, Core::StochasticPropagator< 4 > const &prop_sR,
    Core::StochasticPropagator< 4 > const &prop_dL, Core::StochasticPropagator< 4 > const &prop_sL,
    std::vector< std::pair< Base::Operator, Base::Operator > > const &operator_combinations)
  {

 //   assert(psi1.T() == psi2.T() && psi1.L() == psi2.L());

    std::vector< Core::Correlator< Dirac::Matrix > > threepoints;

    Core::StochasticPropagator< 4 > prop_dR_revert(prop_dR);
    prop_dR_revert.revert();

    Core::StochasticPropagator< 4 > prop_dL_revert(prop_dL);
    prop_dL_revert.revert();

    // loop over all operator combinations
    for (size_t iOp = 0; iOp < operator_combinations.size(); iOp++)
    {

	//RIGHT SIDE     
      Core::StochasticPropagator< 4 > psi1(prop_sR);
      psi1.rightMultiplyOperator(operator_combinations[iOp].first);

      Core::StochasticPropagator< 4 > psi2(prop_dR_revert);
      psi2.rightMultiplyOperator(operator_combinations[iOp].second);

      Core::StochasticPropagator< 4 > O_R(psi1);
      O_R.rightMultiply(psi2);

	//LEFT SIDE

      Core::StochasticPropagator< 4 > psi3(prop_sL);
      psi3.rightMultiplyOperator(operator_combinations[iOp].first);

      Core::StochasticPropagator< 4 > psi4(prop_dL_revert);
      psi4.rightMultiplyOperator(operator_combinations[iOp].second);

      Core::StochasticPropagator< 4 >  O_L(psi3);
      O_L.rightMultiply(psi4);

	// Multiply and take the trace 	

      Core::Correlator< Dirac::Matrix >threepoint(O_R*O_L);

      threepoint.sumOverSpatialVolume(); // this is the projection to zero momentum
      threepoints.push_back(threepoint);
    }

    return threepoints;
  }

}
