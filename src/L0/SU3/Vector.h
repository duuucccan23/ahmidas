#pragma once

#include <L0/SU3/GeneralVector.h>

namespace SU3
{
  template <>
  class GeneralVector< 1 >
  {
    std::complex< double > d_data[3];

    // Some useful constant matrices
    static const GeneralVector< 1 > s_zero;

    public:
      GeneralVector();
      GeneralVector(double const *data);
      GeneralVector(std::complex< double > const *data);
      GeneralVector(GeneralVector< 1 > const &other);
      GeneralVector &operator=(GeneralVector< 1 > const &other);

      void setToRandom();
      void setToZero();

      static GeneralVector< 1 >        random();
      static GeneralVector< 1 > const &zero();

      ~GeneralVector();

      friend std::ostream &operator<<(std::ostream &out, GeneralVector< 1 > const &vec);

      template< typename T >
      GeneralVector< 1 > &operator+=(T const &rhand);

      template< typename T >
      GeneralVector< 1 > &operator-=(T const &rhand);

      template< typename T >
      GeneralVector< 1 > &operator*=(T const &rhand);

      template< typename T >
      GeneralVector< 1 > &operator/=(T const &rhand);

      void leftMultiply(Matrix const &mat);
      void leftMultiply(hcMatrix const &mat);

      std::complex< double > &operator[](short component);
      std::complex< double > const &operator[](short component) const;

      size_t size() const;
  };
  
  typedef class GeneralVector< 1 > Vector;
  
  std::complex< double > innerProduct(Vector const &left, Vector const &right);
}

#include "Vector/Vector.inlines"
