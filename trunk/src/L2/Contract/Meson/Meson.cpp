#include "Meson.ih"

namespace Contract
{

  // this works using the one-end trick and gives all 16 gamma-combinations
  std::vector< Core::Correlator< Dirac::Matrix > > light_meson_twopoint_stochastic(Core::StochasticPropagator< 4 > const &psi1,
                                                                                   Core::StochasticPropagator< 4 > const &psi2)
  {

    assert(psi1.T() == psi2.T() && psi1.L() == psi2.L());

    std::vector< Core::Correlator< Dirac::Matrix > > twopoints;

    Core::StochasticPropagator< 4 > psi2_dagger(psi2);
    psi2_dagger.dagger();

    Dirac::Gamma< 5 >  gamma5;
    Dirac::Gamma< 4 >  gamma0;
    Dirac::Gamma< 1 >  gamma1;
    Dirac::Gamma< 2 >  gamma2;
    Dirac::Gamma< 3 >  gamma3;
    Dirac::Gamma< 45 > gamma0gamma5;
    Dirac::Gamma< 15 > gamma1gamma5;
    Dirac::Gamma< 25 > gamma2gamma5;
    Dirac::Gamma< 35 > gamma3gamma5;
    Dirac::Gamma< 54 > gamma5gamma0;
//     Dirac::Gamma< 51 > gamma5gamma1;
//     Dirac::Gamma< 52 > gamma5gamma2;
//     Dirac::Gamma< 53 > gamma5gamma3;


    Core::StochasticPropagator< 4 > *gamma5_Gamma_psi1;
    Core::StochasticPropagator< 4 > *gamma5_Gamma_psi2_dagger;
    Core::Correlator< Dirac::Matrix > *twopoint;

    // SD: I consider it to be faster to multiply the Gamma to the source indices (Gamma from the right),
    // moreover, this enables us to use the faster '*=' operator instead of '*'

    // note Gamma -> Gamma*gamma5, Gamma' -> gamma5*Gamma'

    // 1) ------ Gamma = gamma5 = Gamma' ------
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    // 2) ------ Gamma = gamma5, Gamma' = gamma0gamma5 ------
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma0;
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    // 3) ------ Gamma = gamma0gamma5, Gamma' = gamma5 ------
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma0;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    // 4) ------ Gamma = gamma0gamma5, Gamma' = gamma0gamma5 ------
    gamma5_Gamma_psi1  = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma0;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma0;
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    // 5) ------ Gamma = gamma0, Gamma' = gamma0 ------
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma0gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma5gamma0;
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    // 6) ------ Gamma = gamma5, Gamma' = gamma0 ------
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma5gamma0;
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    // 7) ------ Gamma = gamma0, Gamma' = gamma5 ------
    gamma5_Gamma_psi1  = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma0gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    // 8) ------ Gamma = gamma0gamma5, Gamma' = gamma0 ------
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma0;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma5gamma0;
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    // 9) ------ Gamma = gamma0, Gamma' = gamma0gamma5 ------
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma0gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma0;
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    // 10) ------ Gamma = gammaigamma0, Gamma' = gammaigamma0 ------ (sum i=1,2,3)
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma1;
    (*gamma5_Gamma_psi1) *= gamma0gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma1;
    (*gamma5_Gamma_psi2_dagger) *= gamma0gamma5;
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma2;
    (*gamma5_Gamma_psi1) *= gamma0gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma2;
    (*gamma5_Gamma_psi2_dagger) *= gamma0gamma5;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma3;
    (*gamma5_Gamma_psi1) *= gamma0gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma3;
    (*gamma5_Gamma_psi2_dagger) *= gamma0gamma5;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;


    // 11) ------ Gamma = gammai, Gamma' = gammai ------ (sum i=1,2,3)
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma1gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma1gamma5;
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma2gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma2gamma5;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma3gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma3gamma5;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    // 12) ------ Gamma = gammaigamma5, Gamma' = gammaigamma5 ------ (sum i=1,2,3)
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma1;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma1;
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma2;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma2;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma3;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma3;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    // 13) ------ Gamma = gammaigamma0, Gamma' = gammai ------ (sum i=1,2,3)
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma1;
    (*gamma5_Gamma_psi1) *= gamma0gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma1gamma5;
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma2;
    (*gamma5_Gamma_psi1) *= gamma0gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma2gamma5;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma3;
    (*gamma5_Gamma_psi1) *= gamma0gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma3gamma5;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    // 14) ------ Gamma = gammai, Gamma' = gammaigamma0 ------ (sum i=1,2,3)
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma1gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma1;
    (*gamma5_Gamma_psi2_dagger) *= gamma0gamma5;
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma2gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma2;
    (*gamma5_Gamma_psi2_dagger) *= gamma0gamma5;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma3gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma3;
    (*gamma5_Gamma_psi2_dagger) *= gamma0gamma5;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    // 15) ------ Gamma = gammaigamma0, Gamma' = gammaigamma5 ------ (sum i=1,2,3)
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma1;
    (*gamma5_Gamma_psi1) *= gamma0gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma1;
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma2;
    (*gamma5_Gamma_psi1) *= gamma0gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma2;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma3;
    (*gamma5_Gamma_psi1) *= gamma0gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma3;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    // 16) ------ Gamma = gammaigamma5, Gamma' = gammaigamma0 ------ (sum i=1,2,3)
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma1;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma1;
    (*gamma5_Gamma_psi2_dagger) *= gamma0gamma5;
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma2;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma2;
    (*gamma5_Gamma_psi2_dagger) *= gamma0gamma5;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma3;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma3;
    (*gamma5_Gamma_psi2_dagger) *= gamma0gamma5;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    // 17) ------ Gamma = gammai, Gamma' = gammaigamma5 ------ (sum i=1,2,3)
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma1gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma1;
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma2gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma2;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma3gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma3;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    // 18) ------ Gamma = gammaigamma5, Gamma' = gammai ------ (sum i=1,2,3)
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma1;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma1gamma5;
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma2;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma2gamma5;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma3;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma3gamma5;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;    

    // 19) ------ Gamma = 1, Gamma' = 1 ------
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma5;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma5;
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    // 20) ------ Gamma = gammaigamma0gamma5, Gamma' = gammaigamma0gamma5 ------ (sum i=1,2,3)
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma1;
    (*gamma5_Gamma_psi1) *= gamma0;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma1;
    (*gamma5_Gamma_psi2_dagger) *= gamma0;
    twopoint = new Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma2;
    (*gamma5_Gamma_psi1) *= gamma0;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma2;
    (*gamma5_Gamma_psi2_dagger) *= gamma0;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    gamma5_Gamma_psi1 = new Core::StochasticPropagator< 4 >(psi1);
    (*gamma5_Gamma_psi1) *= gamma3;
    (*gamma5_Gamma_psi1) *= gamma0;
    gamma5_Gamma_psi2_dagger = new Core::StochasticPropagator< 4 >(psi2_dagger);
    (*gamma5_Gamma_psi2_dagger) *= gamma3;
    (*gamma5_Gamma_psi2_dagger) *= gamma0;
    *twopoint += Core::Correlator< Dirac::Matrix >((*gamma5_Gamma_psi1)*(*gamma5_Gamma_psi2_dagger));
    delete gamma5_Gamma_psi1;
    delete gamma5_Gamma_psi2_dagger;
    twopoint->sumOverSpatialVolume();
    twopoints.push_back(*twopoint);
    delete twopoint;

    return twopoints;
  }

}
