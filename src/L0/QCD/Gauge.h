#pragma once

#include <algorithm>
#include <functional>

#include <L0/SU3/Matrix.h>

namespace QCD
{
  class Gauge
  {
    SU3::Matrix d_data[4];

    public:
      Gauge();
      Gauge(SU3::Matrix const &value);
      Gauge(SU3::Matrix const *values);
      Gauge(Gauge const &other);
      Gauge(double const *data);
      Gauge(std::complex< double > const *data);

      void reunitarize();
      void setToRandom();
      void setToZero();
      void setToIdentity();

      void leftMultiply(Gauge const &other);
      void leftMultiply(SU3::Matrix const &other);

      void rightMultiply(Gauge const &other);
      void rightMultiply(SU3::Matrix const &other);

      bool equals(Gauge const &other, double const relativePrecision) const;

      bool operator==(Gauge const &other) const;

      template< typename T >
      Gauge &operator+=(T const &rhand);

      template< typename T >
      Gauge &operator-=(T const &rhand);

      template< typename T >
      Gauge &operator*=(T const &rhand);

      template< typename T >
      Gauge &operator/=(T const &rhand);

      SU3::Matrix &operator[](short component);
      SU3::Matrix const &operator[](short component) const;

      size_t size() const;
  };

  std::ostream &operator<<(std::ostream &out, Gauge const &gauge);
}

#include "Gauge/Gauge.inlines"
