#ifndef GUARD_QCD_GAUGE_H
#define GUARD_QCD_GAUGE_H

#include <algorithm>
#include <functional>

#include "../../SU3/Matrix/Matrix.h"

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

      void leftMultiply(Gauge const &other);
      void leftMultiply(SU3::Matrix const &other);      

      void rightMultiply(Gauge const &other);      
      void rightMultiply(SU3::Matrix const &other);

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
  };
}

#include "Gauge.inlines"

#endif
