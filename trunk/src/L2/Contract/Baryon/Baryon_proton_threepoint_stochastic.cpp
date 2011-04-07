#include "Baryon.ih"

namespace Contract
{

  std::vector< Core::BaryonCorrelator > proton_threepoint_stochastic(Core::Propagator const &u,
                                                Core::Propagator const &d,
                                                Core::StochasticPropagator <12> const &phi_u,
                                                Core::StochasticPropagator <12> const &phi_d,
                                                Core::StochasticSource <12> const &xi,
                                                size_t t_source, size_t t_sink,
                                                std::vector< Base::Operator > const &my_operators,
                                                Base::BaryonPropagatorProjector const my_projector)
  {

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
    std::vector< Core::Field< Dirac::Matrix > > threepoint
      = u.construct_baryon_with_operator_insertion(d, u, phi_u, phi_d, phi_u, xi,
                                                   Base::bar_PROTON, my_operators,
                                                   t_source, t_sink);

    assert(threepoint.size() == my_operators.size()*2);

    std::vector< Core::BaryonCorrelator > allthreepoints;

    for (size_t opIdx=0; opIdx<my_operators.size(); opIdx++)
    {
      Core::BaryonCorrelator threepoint_UU(threepoint[2*opIdx  ]);
      Core::BaryonCorrelator threepoint_DD(threepoint[2*opIdx+1]);
      threepoint_UU.sumOverSpatialVolume();
      threepoint_DD.sumOverSpatialVolume();

//       threepoint_UU *= my_projector;
//       threepoint_DD *= my_projector;
      allthreepoints.push_back(threepoint_UU);
      allthreepoints.push_back(threepoint_DD);
    }

    return allthreepoints;
  }

  std::vector< Core::BaryonCorrelator > proton_threepoint_stochastic_non_local(Core::Propagator const &u,
                                                Core::Propagator const &d,
                                                Core::StochasticPropagator <12> const &phi_u,
                                                Core::StochasticPropagator <12> const &phi_d,
                                                Core::StochasticSource <12> const &xi,
                                                Core::Field<QCD::Gauge> const &gauge_field,
                                                size_t t_source, size_t t_sink,
                                                std::vector< Base::Operator > const &my_operators,
                                                Base::BaryonPropagatorProjector const my_projector)
  {

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
    std::vector< Core::Field< Dirac::Matrix > > threepoint
      = u.construct_baryon_with_non_local_operator_insertion(d, u, phi_u, phi_d, phi_u, xi, gauge_field,
                                                   Base::bar_PROTON, my_operators,
                                                   t_source, t_sink);

    assert(threepoint.size() == my_operators.size()*2);

    std::vector< Core::BaryonCorrelator > allthreepoints;

    for (size_t opIdx=0; opIdx<my_operators.size(); opIdx++)
    {
      Core::BaryonCorrelator threepoint_UU(threepoint[2*opIdx  ]);
      Core::BaryonCorrelator threepoint_DD(threepoint[2*opIdx+1]);
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
