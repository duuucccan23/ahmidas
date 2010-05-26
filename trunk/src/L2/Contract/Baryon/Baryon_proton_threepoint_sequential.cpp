#include "Baryon.ih"

namespace Contract
{

  std::vector< Core::Correlator > proton_threepoint_sequential(
    Core::Propagator const &bw_prop_u, Core::Propagator const &fw_prop_u,
    Core::Propagator const &bw_prop_d, Core::Propagator const &fw_prop_d,
    Core::Field< QCD::Gauge > * const gauge_field,
    std::vector< Base::Operator > ops)
  {

    size_t const L = bw_prop_u.L();
    size_t const T = bw_prop_u.T();

    assert(L == fw_prop_u.L() && T == fw_prop_u.T());
    assert(L == bw_prop_d.L() && T == bw_prop_d.T());
    assert(L == fw_prop_d.L() && T == fw_prop_d.T());

    std::vector< Core::Correlator > p3p;
    Dirac::Gamma< 5 > gamma5;
//     Dirac::Gamma< 4 > gamma0;

    for (size_t iOp=0; iOp<ops.size(); iOp++)
    {

    {
      Core::Propagator fw_tmp_u(fw_prop_u);
      Core::Propagator bw_tmp_u(bw_prop_u);
      bw_tmp_u.dagger();
      bw_tmp_u *= gamma5;


      bw_tmp_u.multiplyOperator(ops[iOp], gauge_field);

//       switch (ops[iOp])
//       {
//         case Base::op_GAMMA_4:
//           fw_tmp_u.rightMultiply(gamma0);
//           break;
//         case Base::op_UNITY:
//           // nothing to do
//           // fw_tmp_u.rightMultiply(gamma5);
//           break;
//         default:
//         std::cerr << "Error in "
//                   << "std::vector< Core::Correlator > proton_threepoint_sequential(...):\n"
//                   << "Operator with index " << ops[iOp] << " not implemented yet!" << std::endl;
//       }

      Core::Correlator p3p_tmp_u(L, T, bw_tmp_u.contract(fw_tmp_u));
      p3p_tmp_u.sumOverSpatialVolume();
      p3p.push_back(p3p_tmp_u);
    }

    {
      Core::Propagator fw_tmp_d(fw_prop_d);
      Core::Propagator bw_tmp_d(bw_prop_d);
      bw_tmp_d.dagger();
      bw_tmp_d *= gamma5;

      bw_tmp_d.multiplyOperator(ops[iOp], gauge_field);
/*      
//       // since we have to take the transpose later, we replace the above expressions by:
//       bw_tmp_d.rightMultiply(gamma5);
//       bw_tmp_d.conjugate();

      switch (ops[iOp])
      {
        case Base::op_GAMMA_4:
          fw_tmp_d.rightMultiply(gamma0);
          break;
        case Base::op_UNITY:
          // nothing to do
          // fw_tmp_d.rightMultiply(gamma5);
          // fw_tmp_d *= std::complex< double> (-1, 0);
          break;
      }*/

      Core::Correlator p3p_tmp_d(L, T, bw_tmp_d.contract(fw_tmp_d));
      p3p_tmp_d.sumOverSpatialVolume();
      p3p.push_back(p3p_tmp_d);
    }
    }

    return p3p;
  }
}
