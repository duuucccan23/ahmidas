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
        return construct_proton_with_operator_insertion( *this, no2, phi_no1, phi_no2, xi, ops, t_src, t_snk);
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

    Field< QCD::reducedTensor > *field_1;
    Field< QCD::reducedTensor > *field_2;

    for (size_t opNo=0; opNo<ops.size(); opNo++)
    {
  /* under construction */
//      Field< QCD::reducedTensor > *field_1 = new Field< QCD::reducedTensor >(L, T);
//      Field< QCD::reducedTensor > *field_2 = new Field< QCD::reducedTensor >(L, T);



//       Field< QCD::reducedTensor >::iterator It_1(field_1->begin());
//       Field< QCD::reducedTensor >::iterator It_2(field_2->begin());
// 
//       Propagator::const_iterator It_u(S_u.begin());
//       Propagator::const_iterator It_d(S_d.begin());
//       Propagator::const_iterator It_phi_u(phi_u.begin());
//       Propagator::const_iterator It_phi_d(phi_d.begin());
//       Propagator::const_iterator It_xi(xi.begin());
// 
//       size_t count(0);
// 
//       Dirac::Gamma< 5 > gamma5;
// 
//       while(It_1 != field_1->end())
//       {
// 
//         //for spin & color diluted sources left and right multiplication is the same since they're diagonal
//         (*It_phi_d).conjugate();
//         (*It_phi_u).conjugate();
//         (*It_phi_d) =  QCD::Tensor((gamma5 * (*It_phi_d));
//         (*It_phi_u) =  QCD::Tensor((gamma5 * (*It_phi_u));
//         (*It_xi) *= gamma5;
// 
// 
//         phi_d_y = QCD::Tensor((gamma5 * (*It_phi_d));
//         phi_d_y.conjugate();
//         phi_u_y = QCD::Tensor((gamma5 * (*It_phi_u));
//         phi_u_y.conjugate();
//         xi_snk = (*It_xi)*gamma5;
// 
// 
// 
//         QCD::reducedTensor rrU((*It_u), Base::col_RED,   Base::col_RED);
//         QCD::reducedTensor rgU((*It_u), Base::col_RED,   Base::col_GREEN);
//         QCD::reducedTensor rbU((*It_u), Base::col_RED,   Base::col_BLUE);
//         QCD::reducedTensor grU((*It_u), Base::col_GREEN, Base::col_RED);
//         QCD::reducedTensor ggU((*It_u), Base::col_GREEN, Base::col_GREEN);
//         QCD::reducedTensor gbU((*It_u), Base::col_GREEN, Base::col_BLUE);
//         QCD::reducedTensor brU((*It_u), Base::col_BLUE,  Base::col_RED);
//         QCD::reducedTensor bgU((*It_u), Base::col_BLUE,  Base::col_GREEN);
//         QCD::reducedTensor bbU((*It_u), Base::col_BLUE,  Base::col_BLUE);
// 
//         QCD::reducedTensor rrD((*It_d), Base::col_RED,   Base::col_RED);
//         QCD::reducedTensor rgD((*It_d), Base::col_RED,   Base::col_GREEN);
//         QCD::reducedTensor rbD((*It_d), Base::col_RED,   Base::col_BLUE);
//         QCD::reducedTensor grD((*It_d), Base::col_GREEN, Base::col_RED);
//         QCD::reducedTensor gbD((*It_d), Base::col_GREEN, Base::col_BLUE);
//         QCD::reducedTensor ggD((*It_d), Base::col_GREEN, Base::col_GREEN);
//         QCD::reducedTensor brD((*It_d), Base::col_BLUE,  Base::col_RED);
//         QCD::reducedTensor bbD((*It_d), Base::col_BLUE,  Base::col_BLUE);
//         QCD::reducedTensor bgD((*It_d), Base::col_BLUE,  Base::col_GREEN);
// 
//         QCD::reducedTensor rrPhiU((*It_phi_u), Base::col_RED,   Base::col_RED);
//         QCD::reducedTensor rgPhiU((*It_phi_u), Base::col_RED,   Base::col_GREEN);
//         QCD::reducedTensor rbPhiU((*It_phi_u), Base::col_RED,   Base::col_BLUE);
//         QCD::reducedTensor grPhiU((*It_phi_u), Base::col_GREEN, Base::col_RED);
//         QCD::reducedTensor ggPhiU((*It_phi_u), Base::col_GREEN, Base::col_GREEN);
//         QCD::reducedTensor gbPhiU((*It_phi_u), Base::col_GREEN, Base::col_BLUE);
//         QCD::reducedTensor brPhiU((*It_phi_u), Base::col_BLUE,  Base::col_RED);
//         QCD::reducedTensor bgPhiU((*It_phi_u), Base::col_BLUE,  Base::col_GREEN);
//         QCD::reducedTensor bbPhiU((*It_phi_u), Base::col_BLUE,  Base::col_BLUE);
// 
//         QCD::reducedTensor rrPhiD((*It_phi_d), Base::col_RED,   Base::col_RED);
//         QCD::reducedTensor rgPhiD((*It_phi_d), Base::col_RED,   Base::col_GREEN);
//         QCD::reducedTensor rbPhiD((*It_phi_d), Base::col_RED,   Base::col_BLUE);
//         QCD::reducedTensor grPhiD((*It_phi_d), Base::col_GREEN, Base::col_RED);
//         QCD::reducedTensor ggPhiD((*It_phi_d), Base::col_GREEN, Base::col_GREEN);
//         QCD::reducedTensor gbPhiD((*It_phi_d), Base::col_GREEN, Base::col_BLUE);
//         QCD::reducedTensor brPhiD((*It_phi_d), Base::col_BLUE,  Base::col_RED);
//         QCD::reducedTensor bgPhiD((*It_phi_d), Base::col_BLUE,  Base::col_GREEN);
//         QCD::reducedTensor bbPhiD((*It_phi_d), Base::col_BLUE,  Base::col_BLUE);
// 
//         QCD::reducedTensor rrXi((*It_xi), Base::col_RED,   Base::col_RED);
//         QCD::reducedTensor rgXi((*It_xi), Base::col_RED,   Base::col_GREEN);
//         QCD::reducedTensor rbXi((*It_xi), Base::col_RED,   Base::col_BLUE);
//         QCD::reducedTensor grXi((*It_xi), Base::col_GREEN, Base::col_RED);
//         QCD::reducedTensor ggXi((*It_xi), Base::col_GREEN, Base::col_GREEN);
//         QCD::reducedTensor gbXi((*It_xi), Base::col_GREEN, Base::col_BLUE);
//         QCD::reducedTensor brXi((*It_xi), Base::col_BLUE,  Base::col_RED);
//         QCD::reducedTensor bgXi((*It_xi), Base::col_BLUE,  Base::col_GREEN);
//         QCD::reducedTensor bbXi((*It_xi), Base::col_BLUE,  Base::col_BLUE);
// 
// //         std::cout << "count = " << count++ << std::endl;
// //         std::cout << "Xi =\n" << (*It_xi) << std::endl;
// //         std::cout << "Phi_d =\n" << (*It_phi_d) << std::endl;
// //         std::cout << "rrXi = \n" << rrXi << std::endl;
// //         std::cout << "rgXi = \n" << rgXi << std::endl;
// //         std::cout << "rbXi = \n" << rbXi << std::endl;
// //         std::cout << "grXi = \n" << grXi << std::endl;
// //         std::cout << "ggXi = \n" << ggXi << std::endl;
// 
//         (*It_1) = rrD;
//         (*It_2) = rrU;
//         ++It_u;
//         ++It_d;
//         ++It_phi_u;
//         ++It_phi_d;
//         ++It_xi;
//         ++It_1;
//         ++It_2;
//       }

      QCD::reducedTensor sumSink_1_local(std::complex< double >(0.0, 0.0));
      QCD::reducedTensor sumSink_1_global(std::complex< double >(0.0, 0.0));
      QCD::reducedTensor sumSink_2_local(std::complex< double >(0.0, 0.0));
      QCD::reducedTensor sumSink_2_global(std::complex< double >(0.0, 0.0));

      // this part is only summed over sink timeslice
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
        sumSink_1_local += (*field_1)[localIndex];
        sumSink_2_local += (*field_2)[localIndex];
      }
      }
      }
      weave.allReduce(&sumSink_1_local, &sumSink_1_global);
      weave.allReduce(&sumSink_2_local, &sumSink_2_global);

      Field< QCD::reducedTensor >::iterator It_1_new = field_1->begin();
      Field< QCD::reducedTensor >::iterator It_2_new = field_2->begin();
      while(It_1_new != field_1->end())
      {
        (*It_1_new) *= sumSink_1_global;
        (*It_2_new) *= sumSink_2_global;
        ++It_1_new;
        ++It_2_new;
      }

      fields.push_back(field_1);
      fields.push_back(field_2);
    }
    return fields;
  }

}
