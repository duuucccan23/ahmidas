#include "Baryon.ih"

namespace Contract
{
  // keeps source and sink Dirac indices open and writes the sequential sources into seqSrc,
  // which therefore has to be an array of 16 Propagators
  void create_sequential_source_proton_d(Core:: Propagator * const seqSrc,
                                    Core::Propagator const &u1, Core::Propagator const &u2,
                                    Core::Field < QCD::Gauge > &gauge_field,
                                    Smear::fermionFieldSmearing const smearing,
                                    // smearing: iterations and parameter
                                    size_t const nSmear, double const pSmear,
                                    size_t const t_snk)
  {

    size_t const L(u1.L());
    size_t const T(u1.T());
    assert(L == u2.L()        && T == u2.T());
    assert(L == seqSrc[0].L() && T == seqSrc[0].T());

    for (size_t idx = 0; idx < 16; idx++)
      seqSrc[idx] *= 0.0;

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

          QCD::make_sequential_d(tmp, u1[localIndex], u2[localIndex]);
          for (size_t idx = 0; idx < 16; idx++)
          {
            (tmp[idx]).conjugate();
            (tmp[idx]).right_multiply_proton();
            // we absorb the other interpolating field Gamma structure into the source definition
            // for conventional reasons: the contraction then looks the same for the u and d current
            (tmp[idx]).left_multiply_proton();
            (seqSrc[idx])[localIndex] = tmp[idx];
          }
        }
      }
    }

    switch (smearing)
    {
      case Smear::sm_Jacobi:
        for (size_t idx = 0; idx < 16; idx++)
          (seqSrc[idx]).smearJacobi(pSmear, nSmear, gauge_field);
        break;
     default:
       std::cerr << "unknown smearing in function defined in create_sequential_source_proton_d_d.cpp" << std::endl;
       exit(1);
    }

    for (size_t idx = 0; idx < 16; idx++)
      (seqSrc[idx]).rightMultiply(gamma5);
  }
}
