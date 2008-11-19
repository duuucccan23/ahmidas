#pragma once

#include <algorithm>
#include <complex>
#include <functional>
#include <iterator>

#include <L0/Dirac/Identity.h>
#include <L0/Dirac/Pauli.h>
#include <L0/SU3/Vector.h>

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

#include "WeylSpinor/WeylSpinor.inlines"
#include "WeylSpinor/WeylSpinor.templates"
