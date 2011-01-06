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

  std::vector< Core::Field< Dirac::Matrix > * > construct_proton_with_non_local_operator_insertion(
    Propagator const &S_u, Core::Propagator const &S_d,
    StochasticPropagator< 12 > const &phi_u, StochasticPropagator< 12 > const &phi_d,
    StochasticSource< 12 > const &xi, Field<QCD::Gauge> const& gauge_field,
    std::vector< Base::Operator > const &ops,
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

  std::vector< Core::Field< Dirac::Matrix > * >
             Propagator::construct_baryon_with_non_local_operator_insertion(
             Propagator const &no2, Propagator const &no3,
             StochasticPropagator< 12 > const &phi_no1,
             StochasticPropagator< 12 > const &phi_no2,
             StochasticPropagator< 12 > const &phi_no3,
             StochasticSource< 12 > const &xi,
             Field<QCD::Gauge> const& gauge_field,
             Base::BaryonInterpolatingField const iPol,
             std::vector< Base::Operator > const &ops,
             size_t const t_src, size_t const t_snk) const
  {
    switch (iPol)
    {
      case Base::bar_PROTON:
        return construct_proton_with_non_local_operator_insertion(*this, no2, phi_no1, phi_no2, xi, gauge_field, ops, t_src, t_snk);
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

    Dirac::Gamma< 5 > gamma5;
    Dirac::Gamma< 4 > gamma0;
    Dirac::Gamma< 45 > gamma0gamma5;

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
      // Propagator::const_iterator It_xi(xi.begin());

      // just check that we have an operator that is implemented
      switch (ops[opNo])
      {
        case Base::op_GAMMA_4:
        case Base::op_GAMMA_45:
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

      QCD::Tensor tmp;

      //y-dependent loop
      while(It_dd != field_dd->end())
      {

        //for spin & color diluted sources left and right multiplication is the same since they're diagonal

        QCD::Tensor phi_d_y((gamma5 * (*It_phi_d)).dagger());
        QCD::Tensor phi_u_y((gamma5 * (*It_phi_u)).dagger());

        QCD::Tensor S_d_y(*It_d);
        QCD::Tensor S_u_y(*It_u);

        switch (ops[opNo])
        {
          case Base::op_GAMMA_4:
            tmp = gamma0*S_d_y;
            break;
          case Base::op_GAMMA_45:
            tmp = gamma0gamma5*S_d_y;
        }
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

        switch (ops[opNo])
        {
          case Base::op_GAMMA_4:
            tmp = gamma0*S_u_y;
            break;
          case Base::op_GAMMA_45:
            tmp = gamma0gamma5*S_u_y;
           break;
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


      fields.push_back(field_dd);
      fields.push_back(field_uu);
    }
    return fields;
  }




  std::vector< Core::Field< Dirac::Matrix > * > construct_proton_with_non_local_operator_insertion(
             Propagator const &S_u, Propagator const &S_d,
             StochasticPropagator< 12 > const &phi_u,
             StochasticPropagator< 12 > const &phi_d,
             StochasticSource< 12 > const &xi,
             Field<QCD::Gauge> const& gauge_field,
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

    Dirac::Gamma< 5 > gamma5;
    Dirac::Gamma< 4 > gamma0;
    //Dirac::Gamma< 45 > gamma0gamma5;

    Field< QCD::Gauge > gauge_field_shifted4m(gauge_field);
    gauge_field_shifted4m.shift(Base::idx_T, Base::dir_UP);
    // now the gauge field that comes from x-mu is available at x

    StochasticPropagator< 12 > phi_d_shifted4m(phi_d);
    phi_d_shifted4m.shift(Base::idx_T, Base::dir_UP);

    StochasticPropagator< 12 > phi_d_shifted4p(phi_d);
    phi_d_shifted4p.shift(Base::idx_T, Base::dir_DOWN);

    StochasticPropagator< 12 > phi_u_shifted4m(phi_u);
    phi_u_shifted4m.shift(Base::idx_T, Base::dir_UP);

    StochasticPropagator< 12 > phi_u_shifted4p(phi_u);
    phi_u_shifted4p.shift(Base::idx_T, Base::dir_DOWN);

    Propagator S_u_shifted4m(S_u);
    S_u_shifted4m.shift(Base::idx_T, Base::dir_UP);

    Propagator S_u_shifted4p(S_u);
    S_u_shifted4p.shift(Base::idx_T, Base::dir_DOWN);

    Propagator S_d_shifted4m(S_d);
    S_d_shifted4m.shift(Base::idx_T, Base::dir_UP);

    Propagator S_d_shifted4p(S_d);
    S_d_shifted4p.shift(Base::idx_T, Base::dir_DOWN);


    for (size_t opNo=0; opNo<ops.size(); opNo++)
    {
      field_uu = new Field< Dirac::Matrix >(L, T);
      field_dd = new Field< Dirac::Matrix >(L, T);

      Field< Dirac::Matrix >::iterator It_dd(field_dd->begin());
      Field< Dirac::Matrix >::iterator It_uu(field_uu->begin());

      Field< QCD::Gauge >::const_iterator It_g_m(gauge_field_shifted4m.begin());
      Field< QCD::Gauge >::const_iterator It_g(gauge_field.begin());

      Propagator::const_iterator It_u_m(S_u_shifted4m.begin());
      Propagator::const_iterator It_d_m(S_d_shifted4m.begin());
      Propagator::const_iterator It_u_p(S_u_shifted4p.begin());
      Propagator::const_iterator It_d_p(S_d_shifted4p.begin());
      Propagator::const_iterator It_u(S_u.begin());
      Propagator::const_iterator It_d(S_d.begin());

      // note the flavour swap here!!!
      Propagator::const_iterator It_phi_u_m(phi_d_shifted4m.begin());
      Propagator::const_iterator It_phi_d_m(phi_u_shifted4m.begin());
      Propagator::const_iterator It_phi_u_p(phi_d_shifted4p.begin());
      Propagator::const_iterator It_phi_d_p(phi_u_shifted4p.begin());
      Propagator::const_iterator It_phi_u(phi_d.begin());
      Propagator::const_iterator It_phi_d(phi_u.begin());

      //Propagator::const_iterator It_xi(xi.begin());

      // just check that we have an operator that is implemented
      switch (ops[opNo])
      {
        case Base::op_O44:
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

      QCD::Tensor tmp;

      //y-dependent loop
      while(It_dd != field_dd->end())
      {

        //for spin & color diluted sources left and right multiplication is the same since they're diagonal

        QCD::Tensor phi_d_y_m((gamma5 * (*It_phi_d_m)).dagger());
        QCD::Tensor phi_u_y_m((gamma5 * (*It_phi_u_m)).dagger());
        QCD::Tensor phi_d_y_p((gamma5 * (*It_phi_d_p)).dagger());
        QCD::Tensor phi_u_y_p((gamma5 * (*It_phi_u_p)).dagger());
        QCD::Tensor phi_d_y((gamma5 * (*It_phi_d)).dagger());
        QCD::Tensor phi_u_y((gamma5 * (*It_phi_u)).dagger());


        // have to use other gamma for other operators, but for testing it is ok
        phi_d_y_m *= gamma0;
        phi_u_y_m *= gamma0;
        phi_d_y_p *= gamma0;
        phi_u_y_p *= gamma0;
        phi_d_y *= gamma0;
        phi_u_y *= gamma0;


        QCD::Tensor S_d_y_m(*It_d_m);
        QCD::Tensor S_u_y_m(*It_u_m);
        QCD::Tensor S_d_y_p(*It_d_p);
        QCD::Tensor S_u_y_p(*It_u_p);
        QCD::Tensor S_d_y(*It_d);
        QCD::Tensor S_u_y(*It_u);

        {
          QCD::Tensor tmp1(S_d_y_p);
          // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
          tmp1.rightMultiply((*It_g)[Base::idx_T]);
          tmp = tmp1.rightMultiply(phi_d_y);
        }
        {
          QCD::Tensor tmp1(S_d_y_m);
          // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
          tmp1.rightMultiply(((*It_g_m)[Base::idx_T]).dagger());
          tmp -= tmp1.rightMultiply(phi_d_y);
        }
        {
          QCD::Tensor tmp1(phi_d_y_p);
          // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
          tmp1.leftMultiply(((*It_g)[Base::idx_T]).dagger());
          tmp -= tmp1.leftMultiply(S_d_y);
        }
        {
          QCD::Tensor tmp1(phi_d_y_m);
          // part 4 : multiply this by U_mu from right and and shift it up
          tmp1.leftMultiply((*It_g_m)[Base::idx_T]);
          tmp += tmp1.leftMultiply(S_d_y);
        }

        tmp.left_multiply_proton();
        tmp.transposeFull();
        tmp.leftMultiply(dd_part2);

        tmp.getDiracMatrix((*It_dd),  Base::col_RED,   Base::col_RED);
        tmp.getDiracMatrix(colorDiag, Base::col_GREEN, Base::col_GREEN);
        (*It_dd) += colorDiag;
        tmp.getDiracMatrix(colorDiag, Base::col_BLUE,  Base::col_BLUE);
        (*It_dd) += colorDiag;

        {
          QCD::Tensor tmp1(S_u_y_p);
          // part 1 : shift fromRight down and multiply it by U_mu from left ("Field<QCD::Tensor>::rightMultiply")
          tmp1.rightMultiply((*It_g)[Base::idx_T]);
          tmp = tmp1.rightMultiply(phi_u_y);
        }
        {
          QCD::Tensor tmp1(S_u_y_m);
          // part 2 : multiply fromRight by U_mu^dagger from left and shift it up
          tmp1.rightMultiply(((*It_g_m)[Base::idx_T]).dagger());
          tmp -= tmp1.rightMultiply(phi_u_y);
        }
        {
          QCD::Tensor tmp1(phi_u_y_p);
          // part 3 : shift this down and multiply it by U_mu^dagger from right ("Field<QCD::Tensor>::leftMultiply")
          tmp1.leftMultiply(((*It_g)[Base::idx_T]).dagger());
          tmp -= tmp1.leftMultiply(S_u_y);
        }
        {
          QCD::Tensor tmp1(phi_u_y_m);
          // part 4 : multiply this by U_mu from right and and shift it up
          tmp1.leftMultiply((*It_g_m)[Base::idx_T]);
          tmp += tmp1.leftMultiply(S_u_y);
        }

        tmp.rightMultiply(uu_part2);

        tmp.getDiracMatrix((*It_uu),  Base::col_RED,   Base::col_RED);
        tmp.getDiracMatrix(colorDiag, Base::col_GREEN, Base::col_GREEN);
        (*It_uu) += colorDiag;
        tmp.getDiracMatrix(colorDiag, Base::col_BLUE,  Base::col_BLUE);
        (*It_uu) += colorDiag;

        ++It_g_m;
        ++It_g;
        ++It_u_m;
        ++It_d_m;
        ++It_u_p;
        ++It_d_p;
        ++It_u;
        ++It_d;
        ++It_phi_u_m;
        ++It_phi_d_m;
        ++It_phi_u_p;
        ++It_phi_d_p;
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
