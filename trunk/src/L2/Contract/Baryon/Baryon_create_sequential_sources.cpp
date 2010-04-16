#include "Baryon.ih"

namespace Contract
{
//   // this would do the trace over projSrc * sequential_source * projSnk
//   void create_sequential_source_d(Core:: Propagator &seqSrc, Core::Propagator const &u1, Core::Propagator const &u2,
//                              Base::BaryonInterpolatingField const iPol,
//                              BaryonPropagatorProjector const projSrc, BaryonPropagatorProjector const projSnk);

  void create_sequential_source_proton_d(Core:: Propagator * const seqSrc,
                                    Core::Propagator const &u1, Core::Propagator const &u2, size_t const t_snk)
  {

    assert(u1.L() == u2.L() && u1.T() == u2.T() && seqSrc[0].L() == u1.L() && seqSrc[0].T() == u1.T());

    Base::Weave weave(u1.L(), u1.T());

    Dirac::Gamma< 5 > gamma5;
    QCD::Tensor tmp[16];

    size_t localIndex;
    for(size_t idx_Z = 0; idx_Z < u1.L(); idx_Z++)
    {
      for(size_t idx_Y = 0; idx_Y < u1.L(); idx_Y++)
      {
        for(size_t idx_X = 0; idx_X < u1.L(); idx_X++)
        {
          localIndex = weave.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, t_snk);
          /* globalCoordToLocalIndex returns local volume if local data is not available on this cpu */
          if (localIndex == weave.localVolume())
            continue;

          QCD::make_sequential_d(tmp, u1[localIndex], u2[localIndex]);
          for (size_t iDirac=0; iDirac<16; iDirac++)
          {
            (seqSrc[iDirac])[localIndex] = tmp[iDirac];
            (seqSrc[iDirac])[localIndex].conjugate();
            (seqSrc[iDirac])[localIndex] *= gamma5;
            (seqSrc[iDirac])[localIndex].right_multiply_proton();
          }
        }
      }
    }
  }

  void create_sequential_source_proton_u(Core:: Propagator * const seqSrc,
                                    Core::Propagator const &u1, Core::Propagator const &u2)
  {
    assert(false);
  }
}