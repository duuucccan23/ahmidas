#pragma once

#include <L0/SU3/Matrix.h>

#include <complex>

namespace SU3
{
  class Vector
  {
    std::complex< double > d_data[3];

    // Some useful constant matrices
    static const Vector s_zero;

    public:
      Vector();
      Vector(double *data);
      Vector(std::complex< double > *data);
      Vector(Vector const &other);
      Vector &operator=(Vector const &other);

      void setToRandom();
      void setToZero();

      static Vector        random();
      static Vector const &zero();

      ~Vector();

      friend std::ostream &operator<<(std::ostream &out, Vector const &vec);

      template< typename T >
      Vector &operator+=(T const &rhand);

      template< typename T >
      Vector &operator-=(T const &rhand);

      template< typename T >
      Vector &operator*=(T const &rhand);

      template< typename T >
      Vector &operator/=(T const &rhand);

      void leftMultiply(Matrix const &mat);
      void leftMultiply(hcMatrix const &mat);

      std::complex< double > &operator[](short component);
      std::complex< double > const &operator[](short component) const;

      size_t size() const;
  };

  std::complex< double > innerProduct(Vector const &left, Vector const &right);
}

#include "Vector/Vector.inlines"
