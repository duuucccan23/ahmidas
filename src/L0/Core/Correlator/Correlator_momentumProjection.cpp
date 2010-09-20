#include "Correlator.ih"

namespace Core
{
  // performs a summation over timeslices including non-zero momentum projection
  void Correlator::momentumProjection(int const * const momentum)
  {
    size_t const L(this->L());
    size_t const T(this->T());

    if (s_xRelative == NULL || s_xRelative->L() != L ||  s_xRelative->T() != T)
    {
      std::cerr << "Momentum projection not initialized by call of prepareMomentumProjection(...)" << std::endl;
      exit(1);
    }

    std::complex< double > const I(0, 1);

    isolate();
    size_t localIndex;
    for(size_t idx_T = 0;  idx_T < T; idx_T++)
    {
      d_sumTimeslice[idx_T] = Dirac::Matrix(std::complex< double >(0.0, 0.0));
      for(size_t idx_Z = 0; idx_Z < L; idx_Z++)
      {
        for(size_t idx_Y = 0; idx_Y < L; idx_Y++)
        {
          for(size_t idx_X = 0; idx_X < L; idx_X++)
          {

            localIndex = d_weave->globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, idx_T);
            /* globalCoordToLocalIndex returns local volume if local data is not available on this cpu */
            if (localIndex == d_weave->localVolume())
              continue;

            d_sumTimeslice[idx_T] += (*d_data)[localIndex] *
              exp(I * double(std::inner_product(momentum, momentum + 3, (*s_xRelative)[localIndex], 0)));
          }
        }
      }
    }
    d_weave->sumOverTimeSlices(reinterpret_cast< std::complex< double> const * >(d_sumTimeslice),
                               reinterpret_cast< std::complex< double> * >(d_sumTimeslice_global), 16);
  }
}