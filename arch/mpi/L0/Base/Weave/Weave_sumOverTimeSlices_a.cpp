#include "Weave.ih"

namespace Base
{

  void Weave::sumOverTimeSlices(std::complex< double > const *data_send, std::complex< double > *data_recv) const
  {
    // call MPI::Intracomm::Reduce(...)
    (d_grid.grid()).Reduce(static_cast< const void * >(data_send), static_cast< void * >(data_recv),
                           int(d_T), MPI_DOUBLE_COMPLEX, MPI_SUM, 0);
  }

//   template< >
//   void Weave::sumOverTimeSlices(QCD::reducedTensor const *data_send, QCD::reducedTensor *data_recv)
//   {
//     std::complex< double > const *data_s = reinterpret_cast< std::complex< double > const * >(data_send);
//     std::complex< double > *data_r       = reinterpret_cast< std::complex< double > * >(data_recv);
//     // not sure whether this is going to work, if not, sum element by element for all t
//     (d_grid.grid()).Reduce(static_cast< const void * >(data_s), static_cast< void * >(data_r),
//                            int(16*d_T), MPI_DOUBLE_COMPLEX, MPI_SUM, 0);
//   }

}
