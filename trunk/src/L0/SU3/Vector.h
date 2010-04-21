#pragma once

#include <L0/SU3/Tensor.h>

namespace SU3
{
  template <>
  class Tensor< 1 >
  {
    std::complex< double > d_data[3];

    // Some useful constant matrices
    static const Tensor< 1 > s_zero;

    public:
      Tensor();
      Tensor(double const *data);
      Tensor(std::complex< double > const *data);
      Tensor(Tensor< 1 > const &other);
      Tensor &operator=(Tensor< 1 > const &other);

      void setToRandom();
      void setToZero();

      static Tensor< 1 >        random();
      static Tensor< 1 > const &zero();

      ~Tensor();

      friend std::ostream &operator<<(std::ostream &out, Tensor< 1 > const &vec);

      template< typename T >
      Tensor< 1 > &operator+=(T const &rhand);

      template< typename T >
      Tensor< 1 > &operator-=(T const &rhand);

      template< typename T >
      Tensor< 1 > &operator*=(T const &rhand);

      template< typename T >
      Tensor< 1 > &operator/=(T const &rhand);

      void leftMultiply(Matrix const &mat);
      void leftMultiply(hcMatrix const &mat);

      std::complex< double > &operator[](short component);
      std::complex< double > const &operator[](short component) const;

      size_t size() const;
  };
  
  typedef class Tensor< 1 > Vector;
  
  std::complex< double > innerProduct(Vector const &left, Vector const &right);
}

#include "Vector/Vector.inlines"
