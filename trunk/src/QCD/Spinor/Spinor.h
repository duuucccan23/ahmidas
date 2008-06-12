#ifndef GUARD_QCD_SPINOR_H
#define GUARD_QCD_SPINOR_H

#include <algorithm>
#include <complex>

#include "../../Dirac/Gamma/Gamma.h"
#include "../../Dirac/Sigma/Sigma.h"
#include "../../QCD/WeylSpinor/WeylSpinor.h"
#include "../../SU3/Vector/Vector.h"
#include "../../SU3/Matrix/Matrix.h"

namespace QCD
{
  class Spinor
  {
    SU3::Vector d_data[4];

    public:
      Spinor();
      Spinor(SU3::Vector const &value);
      Spinor(Spinor const &other);
      Spinor(double const *data);
      Spinor(std::complex< double > const *data);

      template< typename T >
      Spinor &operator+=(T const &rhand);

      template< typename T >
      Spinor &operator-=(T const &rhand);

      template< typename T >
      Spinor &operator*=(T const &rhand);

      template< typename T >
      Spinor &operator/=(T const &rhand);

      void leftMultiply(SU3::Matrix const &mat);
      void leftMultiply(SU3::hcMatrix const &mat);

      SU3::Vector &operator[](short component);
      SU3::Vector const &operator[](short component) const;

      template< size_t Index >
      void leftMultiply(Dirac::Gamma< Index > const);

      template< size_t Index >
      void leftMultiply(Dirac::Sigma< Index > const);

      WeylSpinor upper();
      WeylSpinor lower();
  };

  template< size_t IndexSink, size_t IndexSource >
  std::complex< double > contractTwoPoint(Dirac::Gamma< IndexSink >, Spinor sink, Dirac::Gamma< IndexSource >, Spinor source);
}

#include "Spinor.inlines"
#include "Spinor.gamma.inlines"

#endif
