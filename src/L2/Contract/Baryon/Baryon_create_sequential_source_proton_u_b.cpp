#include "Baryon.ih"

namespace Contract
{

  // keeps source and sink Dirac indices open and writes the sequential sources into seqSrc,
  // which therefore has to be an array of 16 Propagators
  void create_sequential_source_proton_u(Core:: Propagator * const seqSrc,
                                         Core::Propagator const &d,
                                         Core::Propagator const &u,
                                         size_t const t_snk)
  {
    size_t const L(u.L());
    size_t const T(u.T());
    assert(L == d.L()      && T == d.T());
    assert(L == seqSrc[0].L() && T == seqSrc[0].T());

    Base::Weave weave(L, T);

    Dirac::Gamma< 5 > gamma5;
    QCD::Tensor tmp[16];

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

          QCD::Tensor d_mod(d[localIndex]);
          d_mod.left_multiply_proton();
          d_mod.right_multiply_proton();
          d_mod.transposeFull();

          QCD::make_sequential_u(tmp, d_mod, u[localIndex]);
          for(size_t idx_D = 0; idx_D < 16; idx_D++)
          {
            (seqSrc[idx_D])[localIndex] = QCD::Tensor((tmp[idx_D]).dagger());
            (seqSrc[idx_D])[localIndex].rightMultiply(gamma5);
          }
        }
      }
    }
  }

}
