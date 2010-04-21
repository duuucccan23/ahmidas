#pragma once

#include <L0/SU3/Matrix.h>

#include <complex>

namespace SU3
{
  template< size_t Dim >
  class Tensor
  {
  };
}
//     std::complex< double > *d_data;
// 
//     // Some useful constant matrices
//     static const GeneralVector< Dim > s_zero;
// 
//     public:
//       GeneralVector();
//       GeneralVector(double const *data);
//       GeneralVector(std::complex< double > const *data);
//       GeneralVector(GeneralVector< Dim > const &other);
//       GeneralVector< Dim > &operator=(GeneralVector< Dim > const &other);
// 
//       void setToRandom();
//       void setToZero();
// 
//       static GeneralVector        random();
//       static GeneralVector const &zero();
// 
//       ~GeneralVector();
// 
//       friend std::ostream &operator<<(std::ostream &out, GeneralVector< Dim > const &vec);
// 
//       template< typename T >
//       GeneralVector< Dim > &operator+=(T const &rhand);
// 
//       template< typename T >
//       GeneralVector< Dim > &operator-=(T const &rhand);
// 
//       template< typename T >
//       GeneralVector< Dim > &operator*=(T const &rhand);
// 
//       template< typename T >
//       GeneralVector< Dim > &operator/=(T const &rhand);
// 
//       void leftMultiply(Matrix const &mat, size_t index = 1);
//       void leftMultiply(hcMatrix const &mat, size_t index = 1);
// 
//       std::complex< double > &operator[](short component);
//       std::complex< double > const &operator[](short component) const;
// 
//       size_t size() const;
//   };
// 
//   template< size_t Dim >
//   std::complex< double > innerProduct(GeneralVector< Dim > const &left, GeneralVector< Dim > const &right);
// }
// 
// #include "GeneralVector/GeneralVector.inlines"
