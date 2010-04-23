#include "Baryon.ih"

namespace Contract
{

  // this would do the trace over projector * sequential_source
  void create_sequential_source_proton_u(Core:: Propagator &seqSrc, Core::Propagator const &u1, Core::Propagator const &u2,
                                         size_t const t_snk, Base::BaryonPropagatorProjector const proj)
  {

    assert(u1.L() == u2.L() && u1.T() == u2.T() && seqSrc.L() == u1.L() && seqSrc.T() == u1.T());

    Base::Weave weave(u1.L(), u1.T());

    Dirac::Gamma< 5 > gamma5;
    QCD::Tensor tmp;

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

          QCD::make_sequential_d(tmp, u1[localIndex], u2[localIndex], proj);
          seqSrc[localIndex] = tmp;
          seqSrc[localIndex].conjugate();
          seqSrc[localIndex] *= gamma5;
          seqSrc[localIndex].right_multiply_proton();
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
