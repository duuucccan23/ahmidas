#ifndef GUARD_CORE_CORRELATOR_H
#define GUARD_CORE_CORRELATOR_H

#include <complex>

namespace Core
{
  template< size_t T >
  class Correlator
  {
    std::complex< double > d_data[T];
    
    public:
      Correlator(std::complex< double > *data, size_t sourceSlice);
  };  
}

#include "Correlator_a.template"

#endif
