#include "Propagator.ih"

namespace Core
{
  void Propagator::changeBoundaryConditions_uniformToFixed(size_t timesliceSource, size_t timesliceBoundary)
  {
    Base::Weave weave(L(), T());
    size_t localIndex;
    std::complex< double > phase;
    long idx_T_new;

    //std::cout << "Changing boundary conditions of propagator: uniform to fixed....";
    assert(timesliceSource < T());
    assert(timesliceBoundary < T());
    isolate();

    long const t_s = timesliceSource - (timesliceSource > timesliceBoundary ? 0 : T());
    std::complex< double > const const_phase(exp(std::complex< double >(0, - M_PI * t_s / T())));

    for(size_t idx_T = 0; idx_T < T(); idx_T++)
    {
      idx_T_new = idx_T - (idx_T > timesliceBoundary ? 0 : T());
      phase = const_phase * exp(std::complex< double >(0, M_PI * idx_T_new / T()));
      for(size_t idx_Z = 0; idx_Z < L(); idx_Z++)
      {
        for(size_t idx_Y = 0; idx_Y < L(); idx_Y++)
        {
          for(size_t idx_X = 0; idx_X < L(); idx_X++)
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
    //std::cout << "Done.\n";
    weave.barrier();
  }
}
