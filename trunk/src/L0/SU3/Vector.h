#ifndef GUARD_SU3_VECTOR_H
#define GUARD_SU3_VECTOR_H

#include <L0/SU3/Matrix.h>

#include <complex>

namespace SU3
{
//  std::ostream &operator<<(std::ostream &out, Vector const &vec);
  
  class Vector
  {
    std::complex< double > d_data[3];
    
    // Some useful constant matrices
    static const Vector s_zero;
    static const Vector s_basis[3];

    public:
      Vector();
      Vector(double *data);
      Vector(std::complex< double > *data);
      Vector(Vector const &other);
      Vector &operator=(Vector const &other);

      void setToRandom();
      void setToBasis(size_t idx);
      void setToZero();
      
      static Vector        random();
      static Vector const &basis(size_t idx);
      static Vector const &zero();
      
      ~Vector();

      std::complex< double > &operator()(short idx);
      std::complex< double > const &operator()(short idx) const;

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
  };

  std::complex< double > innerProduct(Vector const &left, Vector const &right);
}

#include "Vector/Vector.inlines"

#endif
