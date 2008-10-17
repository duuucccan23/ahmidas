#ifndef GUARD_CORE_CORRELATOR_H
#define GUARD_CORE_CORRELATOR_H

#include <algorithm>
#include <complex>

#include <L0/Base/Weave.h>

namespace Core
{
  template< size_t T >
  class Correlator
  {
    size_t                 *d_references;
    std::complex< double > *d_data;
    
    public:
      Correlator();
      Correlator(std::complex< double > const &value);
      Correlator(Correlator const &other);

      ~Correlator();
      
      std::complex< double > &operator[](size_t const idx);
      std::complex< double > const &operator[](size_t const idx) const;

      void sumOverTimeSlices();

    private:
      void destroy();
      void isolate();
  };  
}

#include "Correlator/Correlator_Correlator_a.template"
#include "Correlator/Correlator_Correlator_b.template"

#include "Correlator/Correlator_destroy.template"
#include "Correlator/Correlator_isolate.template"

#include "Correlator/Correlator_sumOverTimeSlices.template"

#endif
