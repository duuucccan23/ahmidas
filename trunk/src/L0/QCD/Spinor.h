#pragma once

#include <algorithm>
#include <complex>
#include <iostream>

#include <L0/Base/Base.h>
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
      Spinor &operator=(Spinor const &other);

      template< typename T >
      Spinor &operator+=(T const &rhand);

      template< typename T >
      Spinor &operator-=(T const &rhand);

      template< typename T >
      Spinor &operator*=(T const &rhand);

      template< typename T >
      Spinor &operator/=(T const &rhand);

      bool equals(Spinor const &other, double const relPrecision) const;

      void leftMultiply(SU3::Matrix const &mat);
      void leftMultiply(SU3::hcMatrix const &mat);

      SU3::Vector &operator[](size_t const index);
      SU3::Vector const &operator[](size_t const index) const;

      std::complex< double > &operator()(Base::DiracIndex const dir, Base::ColourIndex const col);
      std::complex< double > const &operator()(Base::DiracIndex const dir, Base::ColourIndex const col) const;

      void leftMultiply(Dirac::Identity const);

      template< size_t Index >
      void leftMultiply(Dirac::Gamma< Index > const);

      template< size_t Index >
      void leftMultiply(Dirac::Sigma< Index > const);

      WeylSpinor upper();
      WeylSpinor lower();

      void setToRandom();
      void setToZero();

      size_t size() const;
      friend std::ostream &operator<<(std::ostream &out, Spinor const &spinor);
  };
  std::ostream &operator<<(std::ostream &out, Spinor const &spinor);

  //Please not that in the following, the left spinor is complex conjugated automatically in this inner product!
  std::complex< double > innerProduct(Spinor const &left, Spinor const &right);
}

#include "Spinor/Spinor.inlines"
#include "Spinor/Spinor.gamma.inlines"
