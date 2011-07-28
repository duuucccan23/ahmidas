#include "Baryon.ih"

namespace Contract
{

  std::vector< Core::BaryonCorrelator > proton_threepoint_stochastic(
                                                Core::Propagator const &S_u,
                                                Core::Propagator const &S_d,
                                                Core::Propagator const &S_u_local,
                                                Core::Propagator const &S_d_local,
                                                Core::StochasticPropagator <12> const &phi_u,
                                                Core::StochasticPropagator <12> const &phi_d,
                                                Core::StochasticSource <12> const &xi_u,
                                                Core::StochasticSource <12> const &xi_d,
                                                size_t t_source, size_t t_sink,
                                                std::vector< Base::Operator > const &my_operators,
                                                std::vector< Base::BaryonPropagatorProjector > const &my_projectors_u,
                                                std::vector< Base::BaryonPropagatorProjector > const &my_projectors_d)
  {

    assert(S_u.L() == S_d.L() && S_u.T() == S_d.T());

    assert(my_projectors_u.size() == my_operators.size());
    assert(my_projectors_d.size() == my_operators.size());

    std::vector< Core::Field< Dirac::Matrix > > fields;
    if (my_operators.size() == 0)
      return std::vector< Core::BaryonCorrelator >();

    size_t const L(S_u.L());
    size_t const T(S_u.T());

    Base::Weave weave(L, T);

    Dirac::Gamma< 5 > gamma5;

    // position labels for each lattice site
    Core::Field< size_t > timeLabel(L , T);
    for(size_t idx_T = 0; idx_T < T; idx_T++)
    {
      size_t localIndex;
      for(size_t idx_Z = 0; idx_Z < L; idx_Z++)
      {
        for(size_t idx_Y = 0; idx_Y < L; idx_Y++)
        {
          for(size_t idx_X = 0; idx_X < L; idx_X++)
          {
            localIndex = weave.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, idx_T);

            if (localIndex == weave.localVolume())
              continue;

            timeLabel[localIndex] = idx_T;
          }
        }
      }
    }


    for (size_t opNo=0; opNo<my_operators.size(); opNo++)
    {

      // additional iterators for summation over sink timeslice
      Core::Propagator::const_iterator It_fw_u(S_u.begin());
      Core::Propagator::const_iterator It_fw_d(S_d.begin());
      Core::Propagator::const_iterator It_xi_u(xi_d.begin());
      Core::Propagator::const_iterator It_xi_d(xi_u.begin());

      // just check that we have an operator that is implemented
      switch (my_operators[opNo])
      {
        case Base::op_GAMMA_4:
        case Base::op_GAMMA_45:
        case Base::op_GAMMA_15:
        case Base::op_GAMMA_25:
        case Base::op_GAMMA_35:
          break;
        default:
        std::cerr << "Error in "
                  << "Contract::proton_threepoint_stochastic(...):\n"
                  << "Operator with index " << my_operators[opNo] << " not implemented yet!" << std::endl;
      }


      std::complex<double> const I (0, 1);
      std::complex<double> const COMPLEX_ZERO (0, 0);

      // standard constructor fills Tensor with zeros
      QCD::Tensor dd_part2_local = QCD::Tensor();
      QCD::Tensor uu_part2_local = QCD::Tensor();

      // this part is only summed over sink timeslice (here one could insert momentum ...)
      for (Core::Field< size_t >::const_iterator It_t = timeLabel.begin(); It_t != timeLabel.end(); ++It_t)
      {

        if (*It_t != t_sink)
        {
          ++It_fw_u;
          ++It_fw_d;
          ++It_xi_u;
          ++It_xi_d;
          continue;
        }

        // we already account for the flavor change in the 
        QCD::Tensor const xi_u_snk((*It_xi_u)*gamma5);
        QCD::Tensor xi_d_snk((*It_xi_d)*gamma5);
        // actually we leave a transposeFull() here since xi is diagonal

        /*
        assert(xi_d_snk(0)   != COMPLEX_ZERO);
        assert(xi_d_snk(1)   == COMPLEX_ZERO);
        assert(xi_d_snk(2)   == COMPLEX_ZERO);
        assert(xi_d_snk(3)   == COMPLEX_ZERO);
        assert(xi_d_snk(7)   == COMPLEX_ZERO);
        assert(xi_d_snk(11)  == COMPLEX_ZERO);
        assert(xi_d_snk(13)  != COMPLEX_ZERO);
        assert(xi_d_snk(143) != COMPLEX_ZERO);
        assert(xi_u_snk(0)   != COMPLEX_ZERO);
        assert(xi_u_snk(1)   == COMPLEX_ZERO);
        assert(xi_u_snk(2)   == COMPLEX_ZERO);
        assert(xi_u_snk(3)   == COMPLEX_ZERO);
        assert(xi_u_snk(7)   == COMPLEX_ZERO);
        assert(xi_u_snk(11)  == COMPLEX_ZERO);
        assert(xi_u_snk(13)  != COMPLEX_ZERO);
        assert(xi_u_snk(143) != COMPLEX_ZERO);
        */

        QCD::Tensor S_d_xf(*It_fw_d);
        S_d_xf.left_multiply_proton();
        S_d_xf.right_multiply_proton();
        S_d_xf.transposeFull();

        QCD::Tensor const S_u_xf(*It_fw_u);

        QCD::Tensor tmp_d_ti_xi;
        QCD::Tensor tmp_u_ti_xi;

        QCD::make_sequential_d(tmp_d_ti_xi, S_u_xf, S_u_xf, my_projectors_d[opNo]);
        QCD::make_sequential_u(tmp_u_ti_xi, S_d_xf, S_u_xf, my_projectors_u[opNo]);


        xi_d_snk.right_multiply_proton();
        xi_d_snk.transposeFull();
        tmp_d_ti_xi.rightMultiply(xi_d_snk);
        tmp_u_ti_xi.leftMultiply(xi_u_snk);

        dd_part2_local += tmp_d_ti_xi;
        uu_part2_local += tmp_u_ti_xi;


        ++It_fw_u;
        ++It_fw_d;
        ++It_xi_u;
        ++It_xi_d;

      }

      assert(It_fw_u == S_u.end());

      
      std::complex< double > tmp_complex_array_local[144];
      std::complex< double > tmp_complex_array[144];

      for(size_t i=0; i<144; i++)
        tmp_complex_array_local[i] = dd_part2_local(i);

      weave.allReduce(tmp_complex_array_local, tmp_complex_array, 144);
      QCD::Tensor const dd_part2(tmp_complex_array);

      for(size_t i=0; i<144; i++)
        tmp_complex_array_local[i] = uu_part2_local(i);

      weave.allReduce(tmp_complex_array_local, tmp_complex_array, 144);
      QCD::Tensor const uu_part2(tmp_complex_array);


      Dirac::Matrix colorDiag;
      QCD::Tensor tmp;

      Core::Field< Dirac::Matrix >field_uu(L, T);
      Core::Field< Dirac::Matrix >field_dd(L, T);
      Core::Field< Dirac::Matrix >::iterator It_dd(field_dd.begin());
      Core::Field< Dirac::Matrix >::iterator It_uu(field_uu.begin());

      Core::Propagator::const_iterator It_u(S_u_local.begin());
      Core::Propagator::const_iterator It_d(S_d_local.begin());
      Core::Propagator::const_iterator It_phi_u(phi_d.begin());
      Core::Propagator::const_iterator It_phi_d(phi_u.begin());

      //y-dependent loop
      while(It_dd != field_dd.end())
      {

        //for spin & color diluted sources left and right multiplication is the same since they're diagonal

        QCD::Tensor phi_d_y((gamma5 * (*It_phi_d)).dagger());
        QCD::Tensor phi_u_y((gamma5 * (*It_phi_u)).dagger());

        QCD::Tensor S_d_y(*It_d);
        QCD::Tensor S_u_y(*It_u);

        switch (my_operators[opNo])
        {
          case Base::op_GAMMA_4:
          {
            Dirac::Gamma< 4 > gamma0;
            tmp = gamma0*S_d_y;
            break;
          }
          case Base::op_GAMMA_45:
          {
            Dirac::Gamma< 45 > gamma0gamma5;
            tmp = gamma0gamma5*S_d_y;
            break;
          }
          case Base::op_GAMMA_15:
          {
            Dirac::Gamma< 15 > gamma1gamma5;
            tmp = gamma1gamma5*S_d_y;
            break;
          }
          case Base::op_GAMMA_25:
          {
            Dirac::Gamma< 25 > gamma2gamma5;
            tmp = gamma2gamma5*S_d_y;
            break;
          }
          case Base::op_GAMMA_35:
          {
            Dirac::Gamma< 35 > gamma3gamma5;
            tmp = gamma3gamma5*S_d_y;
            break;
          }
        }
        tmp.left_multiply_proton();
        tmp.rightMultiply(phi_d_y);
        tmp.transposeFull();
        tmp.leftMultiply(dd_part2);

        tmp.getDiracMatrix((*It_dd),  Base::col_RED,   Base::col_RED);
        tmp.getDiracMatrix(colorDiag, Base::col_GREEN, Base::col_GREEN);
        (*It_dd) += colorDiag;
        tmp.getDiracMatrix(colorDiag, Base::col_BLUE,  Base::col_BLUE);
        (*It_dd) += colorDiag;

        switch (my_operators[opNo])
        {
          case Base::op_GAMMA_4:
          {
            Dirac::Gamma< 4 > gamma0;
            tmp = gamma0*S_u_y;
            break;
          }
          case Base::op_GAMMA_45:
          {
            Dirac::Gamma< 45 > gamma0gamma5;
            tmp = gamma0gamma5*S_u_y;
            break;
          }
          case Base::op_GAMMA_15:
          {
            Dirac::Gamma< 15 > gamma1gamma5;
            tmp = gamma1gamma5*S_u_y;
            break;
          }
          case Base::op_GAMMA_25:
          {
            Dirac::Gamma< 25 > gamma2gamma5;
            tmp = gamma2gamma5*S_u_y;
            break;
          }
          case Base::op_GAMMA_35:
          {
            Dirac::Gamma< 35 > gamma3gamma5;
            tmp = gamma3gamma5*S_u_y;
            break;
          }
        }
        tmp.rightMultiply(phi_u_y);
        tmp.rightMultiply(uu_part2);

        tmp.getDiracMatrix((*It_uu),  Base::col_RED,   Base::col_RED);
        tmp.getDiracMatrix(colorDiag, Base::col_GREEN, Base::col_GREEN);
        (*It_uu) += colorDiag;
        tmp.getDiracMatrix(colorDiag, Base::col_BLUE,  Base::col_BLUE);
        (*It_uu) += colorDiag;

        ++It_u;
        ++It_d;
        ++It_phi_u;
        ++It_phi_d;
        ++It_dd;
        ++It_uu;
      }


      fields.push_back(field_uu);
      fields.push_back(field_dd);
    }

    assert(fields.size() == my_operators.size()*2);

    std::vector< Core::BaryonCorrelator > allthreepoints;

    for (size_t opIdx=0; opIdx<my_operators.size(); opIdx++)
    {
      Core::BaryonCorrelator threepoint_UU(fields[2*opIdx  ]);
      Core::BaryonCorrelator threepoint_DD(fields[2*opIdx+1]);
      threepoint_UU.sumOverSpatialVolume();
      threepoint_DD.sumOverSpatialVolume();

      allthreepoints.push_back(threepoint_UU);
      allthreepoints.push_back(threepoint_DD);
    }

    return allthreepoints;
  }






// **************************************************************************************
// **************************************************************************************
// **************************************************************************************




/*

  std::vector< Core::BaryonCorrelator > proton_threepoint_stochastic_non_local(Core::Propagator const &u,
                                                Core::Propagator const &d,
                                                Core::StochasticPropagator <12> const &phi_u,
                                                Core::StochasticPropagator <12> const &phi_d,
                                                Core::StochasticSource <12> const &xi_u,
                                                Core::StochasticSource <12> const &xi_d,
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
      = u.construct_baryon_with_non_local_operator_insertion(d, u, phi_u, phi_d, phi_u, xi_u, gauge_field,
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

  */
}
