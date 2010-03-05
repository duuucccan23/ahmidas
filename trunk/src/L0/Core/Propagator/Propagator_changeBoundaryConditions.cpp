#include "Propagator.ih"

namespace Core
{

  void Propagator::changeBoundaryConditions_uniformToFixed(size_t timesliceSource, size_t timesliceBoundary)
  {
    size_t T(d_components->T());
    assert(timesliceSource < T);
    assert(timesliceBoundary < T);

    std::cout << "Changing boundary conditions of propagator ... ";
    std::cout.flush();

    Base::Weave weave(L(), T);
    isolate();
    size_t x4, x1, x2, x3, localIndex;

    long const t_s = timesliceSource > timesliceBoundary ? timesliceSource : timesliceSource-T;
    std::complex< double > const const_phase(exp(std::complex< double >(0, -M_PI*double(t_s)/double(T))));

    for(x4=0; x4<T; x4++)
    {
      long const x4_new = x4 > timesliceBoundary ? x4 : x4-T;
      std::complex< double > const phase(const_phase*exp(std::complex< double >(0, M_PI*double(x4_new)/double(T))));
      for(x3=0; x3<L(); x3++)
      {
      for(x2=0; x2<L(); x2++)
      {
      for(x1=0; x1<L(); x1++)
      {
        localIndex = weave.globalCoordToLocalIndex(x1, x2, x3, x4);
        /* globalCoordToLocalIndex returns local volume if local data is not available on this cpu */
        if (localIndex == weave.localVolume())
          continue;

        (*d_components)[localIndex] *= phase;
      }
      }
      }
    }
    std::cout << "done" << std::endl;
  }

}
