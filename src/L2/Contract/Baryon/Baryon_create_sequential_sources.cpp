#include "Baryon.ih"

namespace Contract
{
//   // this would do the trace over projSrc * sequential_source * projSnk
//   void create_sequential_source_d(Core:: Propagator &seqSrc, Core::Propagator const &u1, Core::Propagator const &u2,
//                              Base::BaryonInterpolatingField const iPol,
//                              BaryonPropagatorProjector const projSrc, BaryonPropagatorProjector const projSnk);

  void create_sequential_source_proton_d(Core:: Propagator * const seqSrc,
                                    Core::Propagator const &u1, Core::Propagator const &u2)
  {

    assert(u1.L() == u2.L() && u1.T() == u2.T() && seqSrc[0].L() == u1.L() && seqSrc[0].T() == u1.T());

    Core::Propagator::iterator I[16] = {
      seqSrc[ 0].begin(), seqSrc[ 1].begin(), seqSrc[ 2].begin(), seqSrc[ 3].begin(),
      seqSrc[ 4].begin(), seqSrc[ 5].begin(), seqSrc[ 6].begin(), seqSrc[ 7].begin(),
      seqSrc[ 8].begin(), seqSrc[ 9].begin(), seqSrc[10].begin(), seqSrc[11].begin(),
      seqSrc[12].begin(), seqSrc[13].begin(), seqSrc[14].begin(), seqSrc[15].begin()};

    Core::Propagator::const_iterator I_u1(u1.begin());
    Core::Propagator::const_iterator I_u2(u2.begin());

    Dirac::Gamma< 5 > gamma5;

    QCD::Tensor tmp[16];

    while (I_u1 != u1.end())
    {
      QCD::make_sequential_d(tmp, *I_u1, *I_u2);
      for (size_t iDirac=0; iDirac<16; iDirac++)
      {
        *(I[iDirac]) = tmp[iDirac];
        *(I[iDirac]) *= gamma5;
        (*(I[iDirac])).right_multiply_proton();
        ++I[iDirac];
      }

      ++I_u1;
      ++I_u2;
    }
  }

  void create_sequential_source_proton_u(Core:: Propagator * const seqSrc,
                                    Core::Propagator const &u1, Core::Propagator const &u2)
  {
    assert(false);
  }
}