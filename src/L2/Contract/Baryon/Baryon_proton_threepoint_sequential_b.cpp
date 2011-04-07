#include "Baryon.ih"
#include <L0/Print.h>

namespace Contract
{

  std::vector< Core::BaryonCorrelator > proton_threepoint_sequential(
    Core::Propagator const &bw_prop_u, Core::Propagator const &fw_prop_u,
    Core::Propagator const &bw_prop_d, Core::Propagator const &fw_prop_d,
    Core::Field< QCD::Gauge > * const gauge_field,
    std::string identifier,  std::vector< Core::BaryonCorrelator > * addCorrs)
  {

    bool const twisted_basis(true);

    size_t const L = bw_prop_u.L();
    size_t const T = bw_prop_u.T();

    assert(L == fw_prop_u.L() && T == fw_prop_u.T());
    assert(L == bw_prop_d.L() && T == bw_prop_d.T());
    assert(L == fw_prop_d.L() && T == fw_prop_d.T());

    std::vector< Core::BaryonCorrelator > p3p;
    Dirac::Gamma< 5 > gamma5;

    if (identifier == "local")
    {
      assert(addCorrs == NULL);
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
          Core::Propagator bw_tmp_u(bw_prop_u);
          bw_tmp_u.dagger();
          bw_tmp_u *= gamma5;
          bw_tmp_u.leftMultiplyOperator(ops[iOp], twisted_basis); //check again!!!
          Core::BaryonCorrelator p3p_tmp_u(bw_tmp_u.contract(fw_prop_u));
          //p3p_tmp_u.sumOverSpatialVolume();
          p3p.push_back(p3p_tmp_u);
        }
        // d current
        {
          Core::Propagator bw_tmp_d(bw_prop_d);
          bw_tmp_d.dagger();
          bw_tmp_d *= gamma5;
          bw_tmp_d.leftMultiplyOperator(ops[iOp], twisted_basis);
          Core::BaryonCorrelator p3p_tmp_d(bw_tmp_d.contract(fw_prop_d));
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
              p3p_tmp_d *= -1.0;
              break;
          }
          //p3p_tmp_d.sumOverSpatialVolume();
          p3p.push_back(p3p_tmp_d);
        }
      }
    }
    else if(identifier == "noether")//just placeholder at the moment
    {
      assert(addCorrs != NULL);
      std::vector< Base::Operator > ops;
      ops.push_back(Base::op_CONSERVED_GAMMA_4);
      ops.push_back(Base::op_CONSERVED_GAMMA_1);
      ops.push_back(Base::op_CONSERVED_GAMMA_2);
      ops.push_back(Base::op_CONSERVED_GAMMA_3);
      Core::Field< Dirac::Matrix > axialCurrent(L, T);
      for (size_t iOp=0; iOp<ops.size(); iOp++)
      {
        {
          Core::Propagator bw_tmp_u(bw_prop_u);
          bw_tmp_u.dagger();
          bw_tmp_u *= gamma5;
          Core::BaryonCorrelator p3p_tmp_u(bw_tmp_u.contractWithOperatorInsertion(ops[iOp], gauge_field, fw_prop_u, axialCurrent));
          //p3p_tmp_u.sumOverSpatialVolume();
          p3p.push_back(p3p_tmp_u);
          Core::BaryonCorrelator p3p_tmp_u_5(axialCurrent);
          //p3p_tmp_u_5.sumOverSpatialVolume();
          addCorrs->push_back(p3p_tmp_u_5);
          axialCurrent *= 0.0;
        }
        {
          Core::Propagator bw_tmp_d(bw_prop_d);
          bw_tmp_d.dagger();
          bw_tmp_d *= gamma5;
          Core::BaryonCorrelator p3p_tmp_d(bw_tmp_d.contractWithOperatorInsertion(ops[iOp], gauge_field, fw_prop_d, axialCurrent));
          //p3p_tmp_d.sumOverSpatialVolume();
          p3p.push_back(p3p_tmp_d);
          Core::BaryonCorrelator p3p_tmp_d_5(axialCurrent);
          //p3p_tmp_d_5.sumOverSpatialVolume();
          addCorrs->push_back(p3p_tmp_d_5);
          axialCurrent  *= 0.0;
        }
      }
    }
    else if(identifier == "1D")//just placeholder at the moment
    {
      assert(addCorrs != NULL);
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
      Core::Field< Dirac::Matrix > axialCurrent(L, T);
      for (size_t iOp=0; iOp<ops.size(); iOp++)
      {
        {
          Core::Propagator bw_tmp_u(bw_prop_u);
          bw_tmp_u.dagger();
          bw_tmp_u *= gamma5;
          Core::BaryonCorrelator p3p_tmp_u(bw_tmp_u.contractWithOperatorInsertion(ops[iOp], gauge_field, fw_prop_u, axialCurrent));
          //p3p_tmp_u.sumOverSpatialVolume();
          p3p.push_back(p3p_tmp_u);
          Core::BaryonCorrelator p3p_tmp_u_5(axialCurrent);
          //p3p_tmp_u_5.sumOverSpatialVolume();
          addCorrs->push_back(p3p_tmp_u_5);
          axialCurrent *= 0.0;
        }
        {
          Core::Propagator bw_tmp_d(bw_prop_d);
          bw_tmp_d.dagger();
          bw_tmp_d *= gamma5;
          Core::BaryonCorrelator p3p_tmp_d(bw_tmp_d.contractWithOperatorInsertion(ops[iOp], gauge_field, fw_prop_d, axialCurrent));
          //p3p_tmp_d.sumOverSpatialVolume();
          p3p.push_back(p3p_tmp_d);
          Core::BaryonCorrelator p3p_tmp_d_5(axialCurrent);
          //p3p_tmp_d_5.sumOverSpatialVolume();
          addCorrs->push_back(p3p_tmp_d_5);
          axialCurrent *= 0.0;
        }
      }
    }
    else
    {
      std::string errmsg("Error in Contract::proton_threepoint_sequential(...):\n");
      errmsg.append("unknown value for parameter std::string identifier\n");
      Print(errmsg, std::cerr);
      exit(1);
    }

    return p3p;
  }
}
