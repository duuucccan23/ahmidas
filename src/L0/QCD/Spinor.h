#ifndef GUARD_QCD_SPINOR_H
#define GUARD_QCD_SPINOR_H

#include <algorithm>
#include <complex>

#include <L0/Dirac/Gamma.h>
#include <L0/Dirac/Sigma.h>
#include <L0/QCD/WeylSpinor.h>
#include <L0/SU3/Vector.h>
#include <L0/SU3/Matrix.h>

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

      SU3::Vector &operator[](size_t const index);
      SU3::Vector const &operator[](size_t const index) const;

      std::complex< double > &operator()(Base::DiracIndex const dir, Base::ColourIndex const col);
      std::complex< double > const &operator()(Base::DiracIndex const dir, Base::ColourIndex const col) const;

      template< size_t Index >
      void leftMultiply(Dirac::Gamma< Index > const);

      template< size_t Index >
      void leftMultiply(Dirac::Sigma< Index > const);

      WeylSpinor upper();
      WeylSpinor lower();

      void setToRandom();
      void setToZero();
  };
}

#include "Spinor/Spinor.inlines"
#include "Spinor/Spinor.gamma.inlines"

#endif
