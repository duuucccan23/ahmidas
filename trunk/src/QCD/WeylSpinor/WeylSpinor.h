#ifndef GUARD_QCD_WEYLSPINOR_H
#define GUARD_QCD_WEYLSPINOR_H

#include <algorithm>
#include <complex>
#include <functional>
#include <iterator>

#include "../../Dirac/Identity/Identity.h"
#include "../../Dirac/Pauli/Pauli.h"
#include "../../SU3/Vector/Vector.h"

namespace QCD
{
  class WeylSpinor
  {
    typedef SU3::Vector* iterator;
    typedef std::reverse_iterator< iterator > reverse_iterator;

    SU3::Vector *d_data;
    
    public:
      WeylSpinor(SU3::Vector *data);

      iterator begin();
      iterator end();

      reverse_iterator rbegin();
      reverse_iterator rend();

      WeylSpinor &operator*=(std::complex< double > const &rhand);

      void leftMultiply(Dirac::Identity);

      template< size_t Index >
      void leftMultiply(Dirac::Pauli< Index >);

      template< size_t Index >
      void leftMultiply(Dirac::Pauli< Index >, std::complex< double > const &factor);     

      void swap(WeylSpinor &other);

      void swap(std::complex< double > const &myFac, WeylSpinor other, std::complex< double > const &hisFac);
      void swap(std::complex< double > const &myFac, WeylSpinor other, std::complex< double > const &hisFac, Dirac::Identity);      
      
      template< size_t Index >
      void swap(std::complex< double > const &myFac, WeylSpinor other, std::complex< double > const &hisFac, Dirac::Pauli< Index >);  
  };
}

#include "WeylSpinor.inlines"
#include "WeylSpinor.templates"

#endif
