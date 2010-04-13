#include "Correlator.ih"

namespace Core
{
  // there must be some way to do this more efficiently!
//   void Correlator::sumOverSpatialVolume(size_t const *momentum)
//   {
//     if (momentum[0]==0 && momentum[1]==0 && momentum[2]==0)
//     {
//       sumOverSpatialVolume();
//       return;
//     }
//     isolate();
//     std::complex< double > phase;
//     const size_t zero(0);
//     size_t x[3] = {0, 0, 0};
//     size_t &x1 = x[0];
//     size_t &x2 = x[1];
//     size_t &x3 = x[1];
//     size_t x4;
//     for(x4=0; x4<T; x4++)
//     {
//       d_sumTimeslice[x4] = Dirac::Matrix(std::complex< double >(0.0, 0.0));
//       for(x3=0; x3<L; x3++)
//       {
//       for(x2=0; x2<L; x2++)
//       {
//       for(x1=0; x1<L; x1++)
//       {
//         if (d_weave->isLocallyAvailable(x1, x2, x3, x4))
//         {
//           phase = exp(  std::complex< double >(0, 2.*M_PI/double(L))
//                       * double(std::inner_product(x, x+3, momentum, zero)));
//           d_sumTimeslice[x4] += (*d_data)[d_weave->globalCoordToLocalIndex(x1, x2, x3, x4)] * phase;
//         }
//       }
//       }
//       }
//     }
//     d_weave->sumOverTimeSlices(reinterpret_cast< std::complex< double> const * >(d_sumTimeslice),
//                                reinterpret_cast< std::complex< double> * >(d_sumTimeslice_global), d_sumTimeslice[0].size());
//   }

  void Correlator::sumOverSpatialVolume()
  {
    isolate();
    size_t x4, x1, x2, x3, localIndex;
    for(x4=0; x4 < T(); x4++)
    {
      d_sumTimeslice[x4] = Dirac::Matrix(std::complex< double >(0.0, 0.0));
      for(x3=0; x3<L(); x3++)
      {
      for(x2=0; x2<L(); x2++)
      {
      for(x1=0; x1<L(); x1++)
      {

        localIndex = d_weave->globalCoordToLocalIndex(x1, x2, x3, x4);
        /* globalCoordToLocalIndex returns local volume if local data is not available on this cpu */
        if (localIndex == d_weave->localVolume())
          continue;

        d_sumTimeslice[x4] += (*d_data)[localIndex];
//         std::cout << "index mapping: ("
//                   << x1 << ","
//                   << x2 << ","
//                   << x3 << ","
//                   << x4 << ") => "
//                   << localIndex << std::endl;
      }
      }
      }
    }
    d_weave->sumOverTimeSlices(reinterpret_cast< std::complex< double> const * >(d_sumTimeslice),
                               reinterpret_cast< std::complex< double> * >(d_sumTimeslice_global), 16);
  }


}
