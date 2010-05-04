#include "Baryon.ih"

namespace Contract
{

  std::vector< Core::Correlator > proton_threepoint_stochastic(Core::Propagator const &u,
                                                Core::Propagator const &d,
                                                Core::StochasticPropagator <12> const &phi_u,
                                                Core::StochasticPropagator <12> const &phi_d,
                                                Core::StochasticSource <12> const &xi,
                                                size_t t_source, size_t t_sink,
                                                std::vector< Base::Operator > const &my_operators,
                                                Base::BaryonPropagatorProjector const my_projector)
  {

    // still under construction ...
//     assert(false);

    assert(u.L() == d.L() && u.T() == d.T());
    if (my_projector == Base::proj_PARITY_PLUS_TM)
    {
      // nothing
    }
    else
    {
      std::cerr << "Error in Contract::proton_twopoint(...)";
      std::cerr << "using projector with index" << my_projector << std::endl;
      std::cerr << "THIS IS NOT IMPLEMENTED YET!" << std::endl;
      exit(1);
    }

    // order of d and u in Propagator::construct_baryon is important!
    std::vector< Core::Field< Dirac::Matrix > * > threepoint
      = u.construct_baryon_with_operator_insertion(d, u, phi_u, phi_d, phi_u, xi,
                                                   Base::bar_PROTON, my_operators,
                                                   t_source, t_sink);

    assert(threepoint.size() == my_operators.size()*2);

    std::vector< Core::Correlator > allthreepoints;

    for (size_t opIdx=0; opIdx<my_operators.size(); opIdx++)
    {
      Core::Correlator threepoint_UU(u.L(), u.T(), threepoint[2*opIdx  ]);
      Core::Correlator threepoint_DD(u.L(), u.T(), threepoint[2*opIdx+1]);
      threepoint_UU.sumOverSpatialVolume();
      threepoint_DD.sumOverSpatialVolume();

//       threepoint_UU *= my_projector;
//       threepoint_DD *= my_projector;
      allthreepoints.push_back(threepoint_UU);
      allthreepoints.push_back(threepoint_DD);
    }

    return allthreepoints;
  }
}
