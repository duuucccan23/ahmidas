#include "Meson.ih"

namespace Contract
{

  std::vector< Core::Correlator< Dirac::Matrix > > light_meson_twopoint(
    Core::Propagator const &psi1, Core::Propagator const &psi2,
    std::vector< std::pair< Base::Operator, Base::Operator > > const &operator_combinations)
  {

    assert(psi1.T() == psi2.T() && psi1.L() == psi2.L());

    std::vector< Core::Correlator< Dirac::Matrix > > twopoints;

    Core::Propagator psi2_dagger(psi2);
    psi2_dagger.revert();

    // loop over all operator combinations
    for (size_t iOp = 0; iOp < operator_combinations.size(); iOp++)
    {
      Core::Propagator first_operator_times_psi1(psi1);
      first_operator_times_psi1.rightMultiplyOperator(operator_combinations[iOp].first);

      Core::Propagator second_operator_times_psi2_dagger(psi2_dagger);
      second_operator_times_psi2_dagger.rightMultiplyOperator(operator_combinations[iOp].second);

      Core::Correlator< Dirac::Matrix >twopoint(first_operator_times_psi1 * second_operator_times_psi2_dagger);
      twopoint.sumOverSpatialVolume(); // this is the projection to zero momentum
      twopoints.push_back(twopoint);
    }

    return twopoints;
  }

}
