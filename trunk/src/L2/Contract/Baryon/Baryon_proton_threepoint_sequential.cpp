#include "Baryon.ih"

namespace Contract
{
  std::vector< Core::Correlator > proton_threepoint_sequential(
    Core:: Propagator const &bw_prop_u, Core::Propagator const &fw_prop_u,
    Core:: Propagator const &bw_prop_d, Core::Propagator const &fw_prop_d,
    std::vector< Base::Operator > ops, Base::BaryonPropagatorProjector const my_projector)
  {

    /* under construction */
    assert(false);

    size_t const L = bw_prop_u.L();
    size_t const T = bw_prop_u.T();

    assert(L == fw_prop_u.L() && T == fw_prop_u.T());
    assert(L == bw_prop_d.L() && T ==bw_prop_d.T());
    assert(L == fw_prop_d.L() && T == fw_prop_d.T());

    std::vector< Core::Correlator > p3p;

    Core::Propagator::const_iterator I_bw_u = bw_prop_u.begin();

    Core::Propagator::const_iterator I_fw(fw_prop_u.begin());

    Dirac::Gamma< 5 > gamma5;


    Core::Propagator fw_tmp_u(L, T);
    Core::Propagator bw_tmp_u(L, T);

    Core::Propagator::iterator I_fw_tmp_u(fw_tmp_u.begin());
    Core::Propagator::iterator I_bw_tmp_u(bw_tmp_u.begin());

//     while (I_fw != fw_prop.end())
//     {
// 
//       (*I_fw_tmp) = (*I_fw);
//       (*I_fw_tmp).left_multiply_proton();
// 
// 
//       for (size_t iOp=0; iOp<ops.size(); iOp++)
//       {
//         switch (ops[iOp])
//         {
//           case Base::op_GAMMA_4:
//             break;
//           default:
//           std::cerr << "Error in "
//                     << "std::vector< Core::Field< Dirac::Matrix > * > construct_proton_with_operator_insertion(...):\n"
//                     << "Operator with index " << op << " not implemented yet!" << std::endl;
//         }
//       }
// 
//       ++I_bw_u;
// 
//       (*I_bw_tmp) *= gamma5;
//       ++I_fw;
//       ++I_fw_tmp;
//       ++I_bw_tmp;
//     }

//     Core::Correlator p3p_tmp(L, T, bw_tmp_u*fw_tmp_u);
//     p3p_tmp.sumOverSpatialVolume();
// 
//     p3p.push_back(p3p_tmp);
// 
//     p3p[iOp].sumOverSpatialVolume();

    return p3p;
  }
}
