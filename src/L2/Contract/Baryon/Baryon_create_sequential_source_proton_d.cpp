#include "Baryon.ih"

namespace Contract
{
  // this would do the trace over projector * sequential_source
  void create_sequential_source_proton_d(Core:: Propagator &seqSrc,
                                    Core::Propagator const &u1, Core::Propagator const &u2,
                                    size_t const t_snk, Base::BaryonPropagatorProjector const proj)
  {

    size_t const L(u1.L());
    size_t const T(u1.T());
    assert(L == u2.L()     && T == u2.T());
    assert(L == seqSrc.L() && T == seqSrc.T());

    Base::Weave weave(L, T);

    Dirac::Gamma< 5 > gamma5;
    QCD::Tensor tmp;

    size_t localIndex;
    for(size_t idx_Z = 0; idx_Z < L; idx_Z++)
    {
      for(size_t idx_Y = 0; idx_Y < L; idx_Y++)
      {
        for(size_t idx_X = 0; idx_X < L; idx_X++)
        {
          localIndex = weave.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, t_snk);
          /* globalCoordToLocalIndex returns local volume if local data is not available on this cpu */
          if (localIndex == weave.localVolume())
            continue;

          QCD::make_sequential_d(tmp, u1[localIndex], u2[localIndex], proj);
          tmp *= gamma5;
          tmp = QCD::Tensor(tmp.dagger());
          tmp.right_multiply_proton();
          seqSrc[localIndex] = tmp;
        }
      }
    }
  }

}
