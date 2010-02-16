#include "Meson.ih"

namespace Contract
{

  // this works using the one-end trick and gives all 16 gamma-combinations
  std::vector< Core::Correlator > light_meson_twopoint_stochastic(Core::StochasticPropagator< 4 > const &psi1,
                                                                  Core::StochasticPropagator< 4 > const &psi2)
  {

    assert(psi1.T() == psi2.T() && psi1.L() == psi2.L());

    std::vector< Core::Correlator > twopoints;

    Core::StochasticPropagator< 4 > psi2_dagger(psi2);
    psi2_dagger.dagger();

    Dirac::Gamma< 5 >  gamma5;
    Dirac::Gamma< 4 >  gamma0;
    Dirac::Gamma< 1 >  gamma1;
    Dirac::Gamma< 2 >  gamma2;
    Dirac::Gamma< 3 >  gamma3;
    Dirac::Unity       one;
    Dirac::Gamma< 45 > gamma0gamma5;
    Dirac::Gamma< 15 > gamma1gamma5;
    Dirac::Gamma< 25 > gamma2gamma5;
    Dirac::Gamma< 35 > gamma3gamma5;
    Dirac::Gamma< 54 > gamma5gamma0;
    Dirac::Gamma< 51 > gamma5gamma1;
    Dirac::Gamma< 52 > gamma5gamma2;
    Dirac::Gamma< 53 > gamma5gamma3;
    Dirac::Sigma< 41 > sigma01;
    Dirac::Sigma< 42 > sigma02;
    Dirac::Sigma< 43 > sigma03;
    Dirac::Sigma< 12 > sigma12;
    Dirac::Sigma< 13 > sigma13;
    Dirac::Sigma< 23 > sigma23;

    Core::StochasticPropagator< 4 > *gamma5_Gamma_psi1;
    Core::StochasticPropagator< 4 > *gamma5_Gamma_psi2_dagger;
    Core::Correlator *twopoint;

    // SD: I consider it to be faster to multiply the Gamma to the source indices (Gamma from the right),
    // moreover, this enables us to use the faster '*=' operator instead of '*'

    // note Gamma -> Gamma*gamma5, Gamma' -> gamma5*Gamma'

    // 1) ------ Gamma = gamma5 = Gamma' ------
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    twopoint = new Core::Correlator(psi1.L(), psi2.T(), (*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;

    // 2) ------ Gamma = gamma5, Gamma' = gamma0gamma5 ------
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma0;
    twopoint = new Core::Correlator(psi1.L(), psi2.T(), (*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;

    // 3) ------ Gamma = gamma0gamma5, Gamma' = gamma5 ------
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma0;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    twopoint = new Core::Correlator(psi1.L(), psi2.T(), (*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;

    // 4) ------ Gamma = gamma0gamma5, Gamma' = gamma0gamma5 ------
    gamma5_Gamma_psi1  = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma0;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma0;
    twopoint = new Core::Correlator(psi1.L(), psi2.T(), (*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;

    // 5) ------ Gamma = gamma0, Gamma' = gamma0 ------
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma0gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma5gamma0;
    twopoint = new Core::Correlator(psi1.L(), psi2.T(), (*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;

    // 6) ------ Gamma = gamma5, Gamma' = gamma0 ------
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma5gamma0;
    twopoint = new Core::Correlator(psi1.L(), psi2.T(), (*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;

    // 7) ------ Gamma = gamma0, Gamma' = gamma5 ------
    gamma5_Gamma_psi1  = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma0gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    twopoint = new Core::Correlator(psi1.L(), psi2.T(), (*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;

    // 8) ------ Gamma = gamma0gamma5, Gamma' = gamma0 ------
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma0;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma5gamma0;
    twopoint = new Core::Correlator(psi1.L(), psi2.T(), (*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;

    // 9) ------ Gamma = gamma0gamma5, Gamma' = gamma0gamma5 ------
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma0gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma0;
    twopoint = new Core::Correlator(psi1.L(), psi2.T(), (*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;

    return twopoints;
  }

}
