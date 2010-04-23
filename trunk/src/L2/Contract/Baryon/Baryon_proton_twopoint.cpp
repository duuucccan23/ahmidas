#include "Baryon.ih"

namespace Contract
{

  Core::Correlator proton_twopoint(Core::Propagator const &u1, Core::Propagator const &u2, Core::Propagator const &d,
                                   Base::BaryonPropagatorProjector const projector)
  {
    assert(u1.L() == d.L() && u1.T() == d.T() && u2.L() == d.L() && u2.T() == d.T());
    if (projector == Base::proj_PARITY_PLUS_TM)
    {
      // nothing
    }
    else
    {
      std::cerr << "Error in Contract::proton_twopoint(...)";
      std::cerr << "using projector with index" << projector << std::endl;
      std::cerr << "THIS IS NOT IMPLEMENTED YET!" << std::endl;
      exit(1);
    }

    // order of d and u in Propagator::construct_baryon is important!
    Core::Correlator twopoint(u1.L(), u1.T(),  u1.construct_baryon(d, u2, Base::bar_PROTON));
    twopoint.sumOverSpatialVolume();
//     std::cout << "\nFull proton two point function:\n" << std::endl;
//     for (size_t t=0; t<u1.T(); t++)
//     {
//       std::cout << "t = " << t << std::endl;
//       std::cout << twopoint[t] << std::endl;
//     }
    twopoint *= projector;
    return twopoint;
  }
}
