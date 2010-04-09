#include "Propagator.ih"

namespace Core
{
  // personal note (SD):
  // at some point I should add an additional parameter "size_t const * sink_momentum" for an obvious reason

  std::vector< Core::Field< QCD::reducedTensor > * > construct_proton_with_operator_insertion(
    Propagator const &S_u, Core::Propagator const &S_d,
    StochasticPropagator< 12 > const &phi_u, StochasticPropagator< 12 > const &phi_d,
    StochasticSource< 12 > const &xi, std::vector< Base::Operator > const &ops,
    size_t const t_src, size_t const t_snk);


  std::vector< Core::Field< QCD::reducedTensor > * >
           Propagator::construct_baryon_with_operator_insertion(
             Propagator const &no2, Propagator const &no3,
             StochasticPropagator< 12 > const &phi_no1,
             StochasticPropagator< 12 > const &phi_no2,
             StochasticPropagator< 12 > const &phi_no3,
             StochasticSource< 12 > const &xi,
             Base::BaryonInterpolatingField const iPol,
             std::vector< Base::Operator > const &ops,
             size_t const t_src, size_t const t_snk) const
  {
    switch (iPol)
    {
      case Base::bar_PROTON:
        return construct_proton_with_operator_insertion(*this, no2, phi_no1, phi_no2, xi, ops, t_src, t_snk);
        break;
      default:
      std::cerr << "unknown interpolating field in "
                << "std::vector< Core::Field< QCD::reducedTensor > * > Propagator::construct_baryon_with_operator_insertion(...)!"
                << std::endl;
      std::cerr << "Aborting..." << std::endl;
      exit(1);
    }

  }


  std::vector< Core::Field< QCD::reducedTensor > * > construct_proton_with_operator_insertion(
             Propagator const &S_u, Propagator const &S_d,
             StochasticPropagator< 12 > const &phi_u,
             StochasticPropagator< 12 > const &phi_d,
             StochasticSource< 12 > const &xi,
             std::vector< Base::Operator > const &ops,
             size_t const t_src, size_t const t_snk)
  {
    std::vector< Field< QCD::reducedTensor > * > fields;
    if (ops.size() == 0)
      return fields;

    size_t const L(S_u.L());
    size_t const T(S_u.T());

    Base::Weave weave(L, T);

    Field< QCD::reducedTensor > *field_uu;
    Field< QCD::reducedTensor > *field_dd;


    for (size_t opNo=0; opNo<ops.size(); opNo++)
    {
      field_uu = new Field< QCD::reducedTensor >(L, T);
      field_dd = new Field< QCD::reducedTensor >(L, T);


      Field< QCD::reducedTensor >::iterator It_dd(field_dd->begin());

      Propagator::const_iterator It_u(S_u.begin());
      Propagator::const_iterator It_d(S_d.begin());
      Propagator::const_iterator It_phi_u(phi_u.begin());
      Propagator::const_iterator It_phi_d(phi_d.begin());
      Propagator::const_iterator It_xi(xi.begin());

//       size_t count(0);

      Dirac::Gamma< 5 > gamma5;
      Dirac::Gamma< 4 > gamma0;

      switch (ops[opNo])
      {
        case Base::op_GAMMA_4:
          break;
        default:
        std::cerr << "Error in "
                  << "std::vector< Core::Field< QCD::reducedTensor > * > construct_proton_with_operator_insertion(...):\n"
                  << "Operator with index " << ops[opNo] << " not implemented yet!" << std::endl;
      }

     // ***************************************************************
     // ***************************************************************
     // trace part is still missing !!! (so is the u_bar*Op*u) ********
     // ***************************************************************
     // ***************************************************************

      QCD::Tensor dd_part2_local[16];

      // this part is only summed over sink timeslice (here one could insert momentum ...)
      size_t localIndex;
      for(size_t x3=0; x3<L; x3++)
      {
      for(size_t x2=0; x2<L; x2++)
      {
      for(size_t x1=0; x1<L; x1++)
      {
        localIndex = weave.globalCoordToLocalIndex(x1, x2, x3, t_snk);
        /* globalCoordToLocalIndex returns local volume if local data is not available on this cpu */
        if (localIndex == weave.localVolume())
          continue;

        QCD::Tensor xi_u_snk(xi[localIndex]*gamma5);
        QCD::Tensor xi_d_snk(xi_u_snk);

        xi_d_snk.right_multiply_proton();

        QCD::Tensor S_d_xf(S_d[localIndex]);
        QCD::Tensor S_u_xf(S_u[localIndex]);

        QCD::Tensor tmp_d[16];
        QCD::make_sequential(tmp_d, S_u_xf, S_u_xf);

        for (size_t idx=0; idx<16; idx++)
        {
          tmp_d[idx].rightMultiply(xi_d_snk);
          dd_part2_local[idx] += tmp_d[idx];
        }
      }
      }
      }

      QCD::Tensor dd_part2[16];
      weave.allReduce(dd_part2_local, dd_part2, 16);


      //y-dependent loop
      while(It_dd != field_dd->end())
      {

        //for spin & color diluted sources left and right multiplication is the same since they're diagonal

        QCD::Tensor phi_d_y(gamma5 * (*It_phi_u));
        phi_d_y.conjugate();
        QCD::Tensor phi_u_y(gamma5 * (*It_phi_d));
        phi_u_y.conjugate();
        QCD::Tensor xi_snk((*It_xi)*gamma5);

        QCD::Tensor S_d_y(*It_d);
//         QCD::Tensor S_u_y(*It_u);

        QCD::Tensor tmp(gamma0*S_d_y);
        tmp.rightMultiply(phi_d_y);
        tmp.left_multiply_proton(); // right_multiply has been applied to Xi



        (*It_dd)   =  QCD::reducedTensor(tmp, dd_part2[ 0]);
//         (*It_dd)  +=  QCD::reducedTensor(tmp, dd_part2[ 5]);
//         (*It_dd)  +=  QCD::reducedTensor(tmp, dd_part2[10]);
//         (*It_dd)  +=  QCD::reducedTensor(tmp, dd_part2[15]);

        ++It_u;
        ++It_d;
        ++It_phi_u;
        ++It_phi_d;
        ++It_dd;
      }


      fields.push_back(field_dd);
      fields.push_back(field_uu);
    }
    return fields;
  }

}
