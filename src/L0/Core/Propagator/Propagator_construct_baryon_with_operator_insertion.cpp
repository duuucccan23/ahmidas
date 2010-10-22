#include "Propagator.ih"

namespace Core
{
  // personal note (SD):
  // at some point I should add an additional parameter "size_t const * sink_momentum" for an obvious reason

  std::vector< Core::Field< Dirac::Matrix > * > construct_proton_with_operator_insertion(
    Propagator const &S_u, Core::Propagator const &S_d,
    StochasticPropagator< 12 > const &phi_u, StochasticPropagator< 12 > const &phi_d,
    StochasticSource< 12 > const &xi, std::vector< Base::Operator > const &ops,
    size_t const t_src, size_t const t_snk);


  std::vector< Core::Field< Dirac::Matrix > * >
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
                << "std::vector< Core::Field< Dirac::Matrix > * > Propagator::construct_baryon_with_operator_insertion(...)!"
                << std::endl;
      std::cerr << "Aborting..." << std::endl;
      exit(1);
    }

  }


  std::vector< Core::Field< Dirac::Matrix > * > construct_proton_with_operator_insertion(
             Propagator const &S_u, Propagator const &S_d,
             StochasticPropagator< 12 > const &phi_u,
             StochasticPropagator< 12 > const &phi_d,
             StochasticSource< 12 > const &xi,
             std::vector< Base::Operator > const &ops,
             size_t const t_src, size_t const t_snk)
  {
    std::vector< Field< Dirac::Matrix > * > fields;
    if (ops.size() == 0)
      return fields;

    size_t const L(S_u.L());
    size_t const T(S_u.T());

    Base::Weave weave(L, T);

    Field< Dirac::Matrix > *field_uu;
    Field< Dirac::Matrix > *field_dd;


    for (size_t opNo=0; opNo<ops.size(); opNo++)
    {
      field_uu = new Field< Dirac::Matrix >(L, T);
      field_dd = new Field< Dirac::Matrix >(L, T);


      Field< Dirac::Matrix >::iterator It_dd(field_dd->begin());
      Field< Dirac::Matrix >::iterator It_uu(field_uu->begin());

      Propagator::const_iterator It_u(S_u.begin());
      Propagator::const_iterator It_d(S_d.begin());
      Propagator::const_iterator It_phi_u(phi_d.begin());
      Propagator::const_iterator It_phi_d(phi_u.begin());
      Propagator::const_iterator It_xi(xi.begin());

      Dirac::Gamma< 5 > gamma5;
      Dirac::Gamma< 4 > gamma0;

      // just check that we have an operator that is implemented
      switch (ops[opNo])
      {
        case Base::op_GAMMA_4:
          break;
        default:
        std::cerr << "Error in "
                  << "std::vector< Core::Field< Dirac::Matrix > * > construct_proton_with_operator_insertion(...):\n"
                  << "Operator with index " << ops[opNo] << " not implemented yet!" << std::endl;
      }


      std::complex<double> const I (0, 1);
      std::complex<double> const COMPLEX_ZERO (0, 0);


      QCD::Tensor dd_part2_local, uu_part2_local;
      dd_part2_local *= COMPLEX_ZERO;
      uu_part2_local *= COMPLEX_ZERO;

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

        // actually we leave a transposeFull() here since xi is diagonal
        QCD::Tensor const xi_u_snk(xi[localIndex]*gamma5);
        QCD::Tensor xi_d_snk(xi_u_snk);

//         assert(xi_d_snk[0]   != COMPLEX_ZERO);
//         assert(xi_d_snk[1]   == COMPLEX_ZERO);
//         assert(xi_d_snk[2]   == COMPLEX_ZERO);
//         assert(xi_d_snk[3]   == COMPLEX_ZERO);
//         assert(xi_d_snk[7]   == COMPLEX_ZERO);
//         assert(xi_d_snk[11]  == COMPLEX_ZERO);
//         assert(xi_d_snk[13]  != COMPLEX_ZERO);
//         assert(xi_d_snk[143] != COMPLEX_ZERO);

        QCD::Tensor S_d_xf(S_d[localIndex]);
        S_d_xf.left_multiply_proton();
        S_d_xf.right_multiply_proton();
        S_d_xf.transposeFull();

        QCD::Tensor const S_u_xf(S_u[localIndex]);

        QCD::Tensor tmp_d_ti_xi;
        QCD::Tensor tmp_u_ti_xi;

        QCD::make_sequential_d(tmp_d_ti_xi, S_u_xf, S_u_xf, Base::proj_PARITY_PLUS_TM);
        QCD::make_sequential_u(tmp_u_ti_xi, S_d_xf, S_u_xf, Base::proj_PARITY_PLUS_TM);

        xi_d_snk.right_multiply_proton();
//         xi_d_snk *= -1.0;


        //tmp_d_ti_xi.leftMultiply(xi_d_snk);
        xi_d_snk.transposeFull();
        tmp_d_ti_xi.rightMultiply(xi_d_snk);
        tmp_u_ti_xi.leftMultiply(xi_u_snk);

        dd_part2_local += tmp_d_ti_xi;
        uu_part2_local += tmp_u_ti_xi;

      }
      }
      }

      QCD::Tensor dd_part2;
      weave.allReduce(&dd_part2_local, &dd_part2);
      QCD::Tensor uu_part2;
      weave.allReduce(&uu_part2_local, &uu_part2);

      Dirac::Matrix colorDiag;

      //y-dependent loop
      while(It_dd != field_dd->end())
      {

        //for spin & color diluted sources left and right multiplication is the same since they're diagonal

        QCD::Tensor phi_d_y((gamma5 * (*It_phi_d)).dagger());
        QCD::Tensor phi_u_y((gamma5 * (*It_phi_u)).dagger());

        QCD::Tensor S_d_y(*It_d);
        QCD::Tensor S_u_y(*It_u);

        QCD::Tensor tmp(gamma0*S_d_y);
        tmp.left_multiply_proton();
        tmp.rightMultiply(phi_d_y);
        tmp.transposeFull();
        tmp.leftMultiply(dd_part2);
        //tmp.rightMultiply(dd_part2);

        tmp.getDiracMatrix((*It_dd),  Base::col_RED,   Base::col_RED);
        tmp.getDiracMatrix(colorDiag, Base::col_GREEN, Base::col_GREEN);
        (*It_dd) += colorDiag;
        tmp.getDiracMatrix(colorDiag, Base::col_BLUE,  Base::col_BLUE);
        (*It_dd) += colorDiag;

        tmp = gamma0*S_u_y;
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


      fields.push_back(field_dd);
      fields.push_back(field_uu);
    }
    return fields;
  }

}
