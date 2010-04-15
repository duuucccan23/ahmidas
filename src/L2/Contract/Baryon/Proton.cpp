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


  Core::Correlator proton_threepoint_d(Core:: Propagator const * const bw_prop, Core:: Propagator const &fw_prop, Base::Operator op)
  {

    /* under construction */
    assert(false);
    assert(bw_prop[0].L() == fw_prop.L() && bw_prop[0].T() == fw_prop.T());

    Core::Propagator::const_iterator I_bw[16] = {
      bw_prop[ 0].begin(), bw_prop[ 1].begin(), bw_prop[ 2].begin(), bw_prop[ 3].begin(),
      bw_prop[ 4].begin(), bw_prop[ 5].begin(), bw_prop[ 6].begin(), bw_prop[ 7].begin(),
      bw_prop[ 8].begin(), bw_prop[ 9].begin(), bw_prop[10].begin(), bw_prop[11].begin(),
      bw_prop[12].begin(), bw_prop[13].begin(), bw_prop[14].begin(), bw_prop[15].begin()};

    Core::Propagator::const_iterator I_fw(fw_prop.begin());

    Dirac::Gamma< 5 > gamma5;

    //QCD::Tensor tmp[16];

    Core::Propagator fw_tmp(fw_prop.L(), fw_prop.T());
    Core::Propagator bw_tmp(fw_prop.L(), fw_prop.T());

    Core::Propagator::iterator I_fw_tmp(fw_tmp.begin());
    Core::Propagator::iterator I_bw_tmp(bw_tmp.begin());

    while (I_fw != fw_prop.end())
    {

      (*I_fw_tmp) = (*I_fw);
      (*I_fw_tmp).left_multiply_proton();

      switch (op)
      {
        case Base::op_GAMMA_4:
          break;
        default:
        std::cerr << "Error in "
                  << "std::vector< Core::Field< Dirac::Matrix > * > construct_proton_with_operator_insertion(...):\n"
                  << "Operator with index " << op << " not implemented yet!" << std::endl;
      }

      for (size_t iDirac=0; iDirac<16; iDirac++)
      {
        //*(I_bw[iDirac]) = gamma5 * (*I_bw[iDirac]);
        ++I_bw[iDirac];
        // now here  we have to sandwich operator which kills two of the indices.
      }

      (*I_bw_tmp) *= gamma5;
      ++I_fw;
      ++I_fw_tmp;
      ++I_bw_tmp;
    }

    Core::Correlator p3p(fw_prop.L(),fw_prop.T(), bw_tmp*fw_tmp);

    p3p.sumOverSpatialVolume();

    return p3p;
  }

}
