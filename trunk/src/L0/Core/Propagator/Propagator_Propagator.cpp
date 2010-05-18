#include "Propagator.ih"

inline Core::Propagator::Propagator(Core::Propagator const &other, size_t const timeslice)
  : d_references(other.d_references),
    d_components(other.d_components)
{
  *d_references += 1;
  assert(timeslice > 0 && timeslice < T());
  
  std::complex< double > ZERO(0,0);
  isolate();

  Base::Weave weave(L(), T());

  size_t x1, x2, x3, x4, localIndex;
  for(x4=0; x4 < T(); x4++)
  {
    for(x3=0; x3 < L(); x3++)
    {
      for(x2=0; x2 < L(); x2++)
      {
        for(x1=0; x1 < L(); x1++)
        {

          localIndex = weave.globalCoordToLocalIndex(x1, x2, x3, x4);
          /* globalCoordToLocalIndex returns local volume if local data is not available on this cpu */
          if (x4 == timeslice || localIndex == weave.localVolume())
            continue;

          d_components[localIndex] *= ZERO;
        }
      }
    }
  }
}

