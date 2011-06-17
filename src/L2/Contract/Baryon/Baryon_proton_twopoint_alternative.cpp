#include "Baryon.ih"

namespace Contract
{

  Core::BaryonCorrelator proton_twopoint_alternative(Core::Propagator const &u1, Core::Propagator const &u2, Core::Propagator const &d,
                                   Base::BaryonPropagatorProjector const projector)
  {
    assert(u1.L() == d.L() && u1.T() == d.T() && u2.L() == d.L() && u2.T() == d.T());

    switch(projector)
    {
      case Base::proj_PARITY_PLUS_TM:
      case Base::proj_NO_PROJECTOR:
        break;
      default:
        std::cerr << "Error in Contract::proton_twopoint_alternative(...)";
        std::cerr << "using projector with index" << projector << std::endl;
        std::cerr << "THIS IS NOT IMPLEMENTED YET!" << std::endl;
        exit(1);
    }

    Dirac::Gamma< 5 > gamma5;
    // Core::Propagator u1g5(gamma5*u1);
    // u1g5 *= gamma5;
    // Core::BaryonCorrelator twopoint(u1g5.construct_baryon(d, u2, Base::bar_PROTON_VAR));
    // order of d and u in Propagator::construct_baryon is important!
    
    Core::Field< Dirac::Matrix > tmp_twopoint(u1.construct_baryon(d, u2, Base::bar_PROTON_VAR));
  
    // now the complete Field has to be multiplied by gamma5 from left and right 
    for(Core::Field< Dirac::Matrix >::iterator Id(tmp_twopoint.begin()); Id != tmp_twopoint.end(); ++Id)
    {
      (*Id) = gamma5 * (*Id);
      (*Id) *= gamma5;
    }
    
    Core::BaryonCorrelator twopoint(tmp_twopoint);
    
    twopoint.sumOverSpatialVolume();
    twopoint *= projector;
    return twopoint;
  }
}
