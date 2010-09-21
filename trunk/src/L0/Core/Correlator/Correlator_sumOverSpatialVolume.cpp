#include "Correlator.ih"

namespace Core
{
  template< typename Datatype >
  void Correlator< Datatype >::sumOverSpatialVolume()
  {
    isolate();
    size_t x1, x2, x3, x4, localIndex;
    for(x4=0; x4 < T(); x4++)
    {
      d_sumTimeslice[x4] = Dataype(std::complex< double >(0.0, 0.0));
      for(x3=0; x3 < L(); x3++)
      {
        for(x2=0; x2 < L(); x2++)
        {
          for(x1=0; x1 < L(); x1++)
          {

            localIndex = d_weave->globalCoordToLocalIndex(x1, x2, x3, x4);
            /* globalCoordToLocalIndex returns local volume if local data is not available on this cpu */
            if (localIndex == d_weave->localVolume())
              continue;

            d_sumTimeslice[x4] += (*d_data)[localIndex];
          }
        }
      }
    }
    d_weave->sumOverTimeSlices(reinterpret_cast< std::complex< double> const * >(d_sumTimeslice),
                               reinterpret_cast< std::complex< double> * >(d_sumTimeslice_global),
                               sizeof(Dataype)/sizeof(std::complex< double >));
  }
}
