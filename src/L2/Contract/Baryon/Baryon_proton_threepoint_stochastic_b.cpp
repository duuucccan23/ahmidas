#include "Baryon.ih"
#include <L0/Print.h>

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
                                                Core::Field<QCD::Gauge> const * const gauge_field,
                                                size_t const t_source, size_t const t_sink,
                                                std::string const &operator_flag,
                                                std::vector< Core::BaryonCorrelator > &addCorrs)
  {

    addCorrs.clear();

    if (operator_flag != "local" && operator_flag != "noether" && operator_flag != "1D")
    {
      std::string errmsg("Error in Contract::proton_threepoint_stochastic(...):\n");
      errmsg.append("unknown value for parameter std::string identifier\n");
      Print(errmsg, std::cerr);
      exit(1);
    }

    bool const twisted_basis(true);
    assert(S_u.L() == S_d.L() && S_u.T() == S_d.T());

    std::vector< Core::Field< std::complex< double > > > fields;
    std::vector< Core::Field< std::complex< double > > > fields_g5;

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


    Core::Propagator const * PHI_U = NULL;
    Core::Propagator const * PHI_D = NULL;
    Core::Propagator const * S_U   = NULL;
    Core::Propagator const * S_D   = NULL;

    // for memory saving and for having a cleaner code
    // non-const variables/objects are hidden from function's scope
    {

      // iterators for summation over sink timeslice
      Core::Propagator::const_iterator It_fw_u(S_u.begin());
      Core::Propagator::const_iterator It_fw_d(S_d.begin());
      Core::Propagator::const_iterator It_xi_u(xi_d.begin());
      Core::Propagator::const_iterator It_xi_d(xi_u.begin());

      // std::complex<double> const I (0, 1);
      // std::complex<double> const COMPLEX_ZERO (0, 0);

      QCD::Tensor dd_part2_local[16];
      QCD::Tensor uu_part2_local[16];

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

        QCD::Tensor S_d_xf(*It_fw_d);
        S_d_xf.left_multiply_proton();
        S_d_xf.right_multiply_proton();
        S_d_xf.transposeFull();

        QCD::Tensor const S_u_xf(*It_fw_u);

        QCD::Tensor tmp_d_ti_xi[16];
        QCD::Tensor tmp_u_ti_xi[16];

        QCD::make_sequential_d(tmp_d_ti_xi, S_u_xf, S_u_xf);
        QCD::make_sequential_u(tmp_u_ti_xi, S_d_xf, S_u_xf);

        xi_d_snk.right_multiply_proton();
        xi_d_snk.transposeFull();
        for (size_t sidx=0; sidx<16; sidx++)
        {
          tmp_d_ti_xi[sidx].rightMultiply(xi_d_snk);
          tmp_u_ti_xi[sidx].leftMultiply(xi_u_snk);
          dd_part2_local[sidx] += tmp_d_ti_xi;
          uu_part2_local[sidx] += tmp_u_ti_xi;
        }

        ++It_fw_u;
        ++It_fw_d;
        ++It_xi_u;
        ++It_xi_d;
      }

      assert(It_fw_u == S_u.end());

      std::complex< double > tmp_complex_array_local[16*144];
      std::complex< double > tmp_complex_array[16*144];
      std::vector< QCD::Tensor > dd_part2;
      std::vector< QCD::Tensor > uu_part2;
       
      // ------------------------------------------------------
      
      for(size_t i=0; i<16*144; i++)
        tmp_complex_array_local[i] = dd_part2_local[i/144](i%144);

      weave.allReduce(tmp_complex_array_local, tmp_complex_array, 16*144);
      for(size_t i=0; i<16; i++)
        dd_part2.push_back(QCD::Tensor(&(tmp_complex_array[i*144])));

      // ------------------------------------------------------
      
      for(size_t i=0; i<16*144; i++)
        tmp_complex_array_local[i] = uu_part2_local[i/144](i%144);

      weave.allReduce(tmp_complex_array_local, tmp_complex_array, 16*144);
      for(size_t i=0; i<16; i++)
        uu_part2.push_back(QCD::Tensor(&(tmp_complex_array[i*144])));

      // ------------------------------------------------------
      
      std::vector< QCD::Tensor > dd_part2_T(dd_part2);
      for(size_t i=0; i<16; i++)
        dd_part2_T[i].transposeFull();

      Core::Propagator PHI_U_tmp((gamma5 * phi_d).dagger());
      Core::Propagator PHI_D_tmp((gamma5 * phi_u).dagger());

      for (Core::Propagator::iterator It_phi_u=PHI_U_tmp.begin(); It_phi_u != PHI_U_tmp.end(); ++It_phi_u)
        (*It_phi_u).rightMultiply(uu_part2);

      for (Core::Propagator::iterator It_phi_d=PHI_D_tmp.begin(); It_phi_d != PHI_D_tmp.end(); ++It_phi_d)
        (*It_phi_d).rightMultiply(dd_part2_T);

      PHI_U = new Core::Propagator(PHI_U_tmp);
      PHI_D = new Core::Propagator(PHI_D_tmp);

      Core::Propagator S_D_tmp(S_d_local);
      for (Core::Propagator::iterator It_S_d = S_D_tmp.begin(); It_S_d != S_D_tmp.end(); ++It_S_d)
        (*It_S_d).left_multiply_proton();
      S_D = new Core::Propagator(S_D_tmp);
      S_U = new Core::Propagator(S_u_local);
    }

    if (operator_flag == "local")
    {
      std::vector< Base::HermitianBilinearOperator > ops;
      ops.push_back(Base::op_G_0);
      ops.push_back(Base::op_G_1);
      ops.push_back(Base::op_G_2);
      ops.push_back(Base::op_G_3);
      ops.push_back(Base::op_G_4);
      ops.push_back(Base::op_G_5);
      ops.push_back(Base::op_G_6);
      ops.push_back(Base::op_G_7);
      ops.push_back(Base::op_G_8);
      ops.push_back(Base::op_G_9);
      ops.push_back(Base::op_G_10);
      ops.push_back(Base::op_G_11);
      ops.push_back(Base::op_G_12);
      ops.push_back(Base::op_G_13);
      ops.push_back(Base::op_G_14);
      ops.push_back(Base::op_G_15);
      for (size_t iOp=0; iOp<ops.size(); iOp++)
      {
        // u current
        {
          Core::Propagator PHI_U_copy(*PHI_U);
          PHI_U_copy.leftMultiplyOperator(ops[iOp], twisted_basis); //check again!!!
          fields.push_back(PHI_U_copy.contract(*S_U));
        }
        // d current
        {
          Core::Propagator PHI_D_copy(*PHI_D);
          PHI_D_copy.leftMultiplyOperator(ops[iOp], twisted_basis);
          Core::Field< Dirac::Matrix > dCurrent(PHI_D_copy.contract(*S_D));
          // this comes from a factor tau_3 for operators in twisted basis
          switch (ops[iOp])
          {
            case Base::op_G_0:
            case Base::op_G_5:
            case Base::op_G_6:
            case Base::op_G_7:
            case Base::op_G_8:
            case Base::op_G_13:
            case Base::op_G_14:
            case Base::op_G_15:
              dCurrent *= -1.0;
              break;
          }
          fields.push_back(dCurrent);
        }
      }
    }
    else if (operator_flag == "noether")
    {
      std::vector< Base::Operator > ops;
      ops.push_back(Base::op_CONSERVED_GAMMA_4);
      ops.push_back(Base::op_CONSERVED_GAMMA_1);
      ops.push_back(Base::op_CONSERVED_GAMMA_2);
      ops.push_back(Base::op_CONSERVED_GAMMA_3);
      Core::Field< Dirac::Matrix > axialCurrent(L, T);
      for (size_t opNo=0; opNo<ops.size(); opNo++)
      {
        //y-dependent loop is realized on "Field" level
        Core::Propagator PHI_U_copy(*PHI_U);
        Core::Propagator PHI_D_copy(*PHI_D);
        fields.push_back(PHI_U_copy.contractWithOperatorInsertion(ops[opNo], gauge_field, *S_U, axialCurrent));
        fields_g5.push_back(axialCurrent);
        fields.push_back(PHI_D_copy.contractWithOperatorInsertion(ops[opNo], gauge_field, *S_D, axialCurrent));
        fields_g5.push_back(axialCurrent);
      }
    }
    else if(operator_flag == "1D")
    {
      std::vector< Base::Operator > ops;
      ops.push_back(Base::op_O44);
      ops.push_back(Base::op_O41);
      ops.push_back(Base::op_O42);
      ops.push_back(Base::op_O43);
      ops.push_back(Base::op_O14);
      ops.push_back(Base::op_O11);
      ops.push_back(Base::op_O12);
      ops.push_back(Base::op_O13);
      ops.push_back(Base::op_O24);
      ops.push_back(Base::op_O21);
      ops.push_back(Base::op_O22);
      ops.push_back(Base::op_O23);
      ops.push_back(Base::op_O34);
      ops.push_back(Base::op_O31);
      ops.push_back(Base::op_O32);
      ops.push_back(Base::op_O33);
      Core::Field< Dirac::Matrix > polarizedCurrent(L, T);
      for (size_t opNo=0; opNo<ops.size(); opNo++)
      {
        //y-dependent loop is realized on "Field" level
        Core::Propagator PHI_U_copy(*PHI_U);
        Core::Propagator PHI_D_copy(*PHI_D);
        fields.push_back(PHI_U_copy.contractWithOperatorInsertion(ops[opNo], gauge_field, *S_U, polarizedCurrent));
        fields_g5.push_back(polarizedCurrent);
        fields.push_back(PHI_D_copy.contractWithOperatorInsertion(ops[opNo], gauge_field, *S_D, polarizedCurrent));
        fields_g5.push_back(polarizedCurrent);
      }
    }

    assert(fields.size() % 2 == 0);

    std::vector< Core::BaryonCorrelator > allthreepoints;

    for (size_t idx=0; idx<fields.size(); idx+=2)
    {
      Core::BaryonCorrelator threepoint_UU(fields[idx  ]);
      Core::BaryonCorrelator threepoint_DD(fields[idx+1]);
      threepoint_UU.sumOverSpatialVolume();
      threepoint_DD.sumOverSpatialVolume();

      allthreepoints.push_back(threepoint_UU);
      allthreepoints.push_back(threepoint_DD);
    }

    assert(fields_g5.size() % 2 == 0);

    for (size_t idx=0; idx<fields_g5.size(); idx+=2)
    {
      Core::BaryonCorrelator threepoint_UU(fields_g5[idx]);
      Core::BaryonCorrelator threepoint_DD(fields_g5[idx+1]);
      threepoint_UU.sumOverSpatialVolume();
      threepoint_DD.sumOverSpatialVolume();
      addCorrs.push_back(threepoint_UU);
      addCorrs.push_back(threepoint_DD);
    }

    delete PHI_U;
    PHI_U = NULL;
    delete PHI_D;
    PHI_D = NULL;
    delete S_U;
    S_U   = NULL;
    delete S_D;
    S_D   = NULL;

    return allthreepoints;
  }

}

