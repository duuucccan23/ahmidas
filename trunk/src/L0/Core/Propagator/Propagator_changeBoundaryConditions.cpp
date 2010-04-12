#include "Propagator.ih"

namespace Core
{

  void Propagator::changeBoundaryConditions_uniformToFixed(size_t timesliceSource, size_t timesliceBoundary)
  {
    size_t T(d_components->T());
    size_t L(d_components->L());

    assert(timesliceSource < T);
    assert(timesliceBoundary < T);

    std::cout << "Changing boundary conditions of propagator: uniform to fixed.\n";

    Base::Weave weave(L, T);
    isolate();
    size_t idx_T, idx_X, idx_Y, idx_Z, localIndex;

    long const t_s = timesliceSource > (timesliceBoundary ? timesliceSource : timesliceSource - T);
    std::complex< double > const const_phase(exp(std::complex< double >(0, - M_PI * double(t_s) / double(T))));

    for(idx_T=0; idx_T<T; idx_T++)
    {
      long const idx_T_new = idx_T > (timesliceBoundary ? idx_T : idx_T - T);
      std::complex< double > const phase(const_phase*exp(std::complex< double >(0, M_PI*double(idx_T_new)/double(T))));
      for(idx_Z=0; idx_Z < L; idx_Z++)
      {
        for(idx_Y=0; idx_Y < L; idx_Y++)
        {
          for(idx_X=0; idx_X < L; idx_X++)
          {
            localIndex = weave.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, idx_T);
            /* globalCoordToLocalIndex returns local volume if local data is not available on this cpu */
            if (localIndex == weave.localVolume())
              continue;

            (*d_components)[localIndex] *= phase;
          }
        }
      }
    }
    std::cout << "Done.\n";
  }
}
