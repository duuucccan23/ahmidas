#include "Baryon.ih"
#include <L0/Base/Weave.h>
#include <L0/QCD/Tensor.h>

namespace Contract
{
  std::vector< Core::Correlator > proton_threepoint_stochastic_naive(Core::Propagator const &u,
                                                          Core::Propagator const &d,
                                                          Core::StochasticPropagator<12> const &u_stoch_at_sink,
                                                          Core::StochasticPropagator<12> const &d_stoch_at_sink,
                                                          Core::StochasticSource<12> const &xi_at_sink,
                                                          Core::Field < QCD::Gauge > const * const gauge_field,
                                                          std::vector< Base::Operator > const &ops,
                                                          size_t const t_src, size_t const t_snk)
  {
    assert(u.L() == d.L() && u.T() == d.T());
    if(gauge_field != NULL)
      assert(gauge_field->L() == d.L() && gauge_field->T() == d.T());

    size_t const L(d.L());
    size_t const T(d.T());

    std::vector< Core::Correlator > threepoints_all;

    if (ops.size() == 0)
      return threepoints_all;;

    Base::Weave weave(L, T);

    Dirac::Gamma< 5 > gamma5;
    Dirac::Gamma< 4 > gamma0;

    Core::Field< Dirac::Matrix > *field_dd;
    Core::Field< Dirac::Matrix > *field_uu;

    for (size_t opNo=0; opNo<ops.size(); opNo++)
    {
      field_dd = new Core::Field< Dirac::Matrix > (L, T);
      field_uu = new Core::Field< Dirac::Matrix > (L, T);

      Core::Field< Dirac::Matrix >::iterator It_dd(field_dd->begin());
      Core::Field< Dirac::Matrix >::iterator It_uu(field_uu->begin());

      Core::Propagator::const_iterator It_u(u.begin());
      Core::Propagator::const_iterator It_d(d.begin());
      Core::Propagator::const_iterator It_phi_d(u_stoch_at_sink.begin());
      Core::Propagator::const_iterator It_phi_u(d_stoch_at_sink.begin());

      size_t pos_snk[4];
      pos_snk[Base::idx_T] = t_snk;

      Dirac::Matrix dd_tmp, uu_tmp;

      size_t count(0);

      //y-dependent loop
      while(It_dd != field_dd->end())
      {
        // note that *It_dd is initialized to zero by default thanks to the default constructor

        // this part is only summed over sink timeslice (here one could in principle insert momentum ...)
        for(size_t idx_X=0; idx_X<L; idx_X++)
        {
          pos_snk[Base::idx_X] = idx_X;
          for(size_t idx_Y=0; idx_Y<L; idx_Y++)
          {
            pos_snk[Base::idx_Y] = idx_Y;
            for(size_t idx_Z=0; idx_Z<L; idx_Z++)
            {
              pos_snk[Base::idx_Z] = idx_Z;

              QCD::Tensor const xi_snk(xi_at_sink(pos_snk));

              QCD::Tensor const tmp_d_from_source(*It_d);
              QCD::Tensor const tmp_u_from_source(*It_u);
//               if (count == 65 && idx_X+idx_Y+idx_Z == 0)
//               {
//                   std::cout << tmp_d_from_source << std::endl;
//               }
              // note: flavour change introduced by hermiticity trick in twisted mass
              // has already been accounted for in declaration of iterators, see above!

              QCD::Tensor tmp_d_from_sink(*It_phi_d);
              QCD::Tensor tmp_u_from_sink(*It_phi_u);

              tmp_d_from_sink.leftMultiplySpinColorDilutedConj(xi_snk);
              tmp_u_from_sink.leftMultiplySpinColorDilutedConj(xi_snk);


              // note that the following would do the same but more inefficiently,
              // since dilution is not accounted for:
              //  xi_d_snk.conjugate();
              //  xi_u_snk.conjugate();
              //  tmp_d_from_sink.leftMultiply(xi_d_snk);
              //  tmp_u_from_sink.leftMultiply(xi_u_snk);


//               if (idx_X==1 && idx_Y==2 && idx_Z==3)
//               {
//                 std::cout << xi_u_snk << std::endl;
//                 std::cout << *It_phi_d << std::endl;
//                 std::cout << tmp_d_from_sink << std::endl;
//                 exit(1);
//               }

              // important note: QCD::Tensor::dagger() does not change the QCD::Tensor it is applied to but returns a QCD::hcTensor

              // apply gamma5 hermiticity trick
              tmp_d_from_sink = QCD::Tensor(tmp_d_from_sink.dagger());
              tmp_d_from_sink.rightMultiply(gamma5);
              tmp_d_from_sink *= gamma5;
              tmp_u_from_sink = QCD::Tensor(tmp_u_from_sink.dagger());
              tmp_u_from_sink.rightMultiply(gamma5);
              tmp_u_from_sink *= gamma5;


              switch (ops[opNo])
              {
                case Base::op_GAMMA_4:
                  tmp_d_from_sink *= gamma0;
                  tmp_u_from_sink *= gamma0;
                  break;
                case Base::op_UNITY:
                  // nothing  to do
                  // tmp_d_from_sink *= gamma5;
                  // tmp_u_from_sink *= gamma5;
                  // tmp_d_from_sink *= std::complex< double >(-1, 0);
                  break;
                default:
                std::cerr << "Error in "
                          << "std::vector< Core::Field< Dirac::Matrix > * > construct_proton_with_operator_insertion(...):\n"
                          << "Operator with index " << ops[opNo] << " not implemented yet!" << std::endl;
              }



              tmp_d_from_sink.leftMultiply(tmp_d_from_source);
              tmp_u_from_sink.leftMultiply(tmp_u_from_source);

              QCD::Tensor S_d_xf(d(pos_snk));
              QCD::Tensor S_u_xf(u(pos_snk));

              // now the twopoint routine should do the job ...
              getDiracMatrix(dd_tmp, S_u_xf, tmp_d_from_sink, S_u_xf, Base::bar_PROTON);
              (*It_dd) += dd_tmp;
              getDiracMatrix(uu_tmp, S_u_xf, S_d_xf, tmp_u_from_sink, Base::bar_PROTON);
              (*It_uu) += uu_tmp;
              getDiracMatrix(uu_tmp, tmp_u_from_sink, S_d_xf, S_u_xf, Base::bar_PROTON);
              (*It_uu) += uu_tmp;

            }
          }
        }

        ++It_u;
        ++It_d;
        ++It_phi_u;
        ++It_phi_d;
        ++It_dd;
        ++It_uu;
        ++count;
      }

      Core::Correlator tp_dd(u.L(), u.T(), field_dd);
      Core::Correlator tp_uu(u.L(), u.T(), field_uu);
      tp_dd.sumOverSpatialVolume();
      tp_uu.sumOverSpatialVolume();
      threepoints_all.push_back(tp_dd);
      threepoints_all.push_back(tp_uu);
    }

    return threepoints_all;
  }
}
