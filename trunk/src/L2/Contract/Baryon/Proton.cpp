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

  std::vector< Core::Correlator > proton_threepoint_stochastic(Core::Propagator const &u,
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
    std::vector< Core::Field< Dirac::Matrix > * > threepoint
      = u.construct_baryon_with_operator_insertion(d, u, phi_u, phi_d, phi_u, xi,
                                                   Base::bar_PROTON, my_operators,
                                                   t_source, t_sink);

    assert(threepoint.size() == my_operators.size()*2);

    std::vector< Core::Correlator > allthreepoints;

    std::cout << "\nFull proton two point function:\n" << std::endl;
    for (size_t opIdx=0; opIdx<my_operators.size(); opIdx++)
    {
      Core::Correlator threepoint_UU(u.L(), u.T(), threepoint[2*opIdx  ]);
      Core::Correlator threepoint_DD(u.L(), u.T(), threepoint[2*opIdx+1]);
      threepoint_UU.sumOverSpatialVolume();
      threepoint_DD.sumOverSpatialVolume();
      std::cout << "\n u_bar * operator * u  (operator no: " << my_operators[opIdx] << ")\n" << std::endl;
      for (size_t t=0; t<u.T(); t++)
      {
        std::cout << "t = " << t << std::endl;
        std::cout << threepoint_UU[t] << std::endl;
      }
      std::cout << "\n d_bar * operator * d  (operator no: " << my_operators[opIdx] << ")\n" << std::endl;
      for (size_t t=0; t<u.T(); t++)
      {
        std::cout << "t = " << t << std::endl;
        std::cout << threepoint_DD[t] << std::endl;
      }
      std::cout << std::endl;
      threepoint_UU *= my_projector;
      threepoint_DD *= my_projector;
      allthreepoints.push_back(threepoint_UU);
      allthreepoints.push_back(threepoint_DD);
    }

    return allthreepoints;
  }

}
