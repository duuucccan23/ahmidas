#pragma once

#include <complex>
#include <algorithm>
#include <functional>

#include <L0/Base/Base.h>
#include <L0/Dirac/Gamma.h>


namespace Dirac
{
  // outer_product returns essentially 16 Matrices, 
  // but the index order has to be defined
  // (Result[a,4*b])[c,4*d] = ?
  enum OuterProductIndexOrder
  {
    order_OUTER_FIXED,        // A_ac*B_db
    order_FIRST_FIXED,        // A_ab*B_cd
    order_FIRST_OUTER_DELTA,  // delta_ad*[A*B]_cb
    order_SECOND_OUTER_DELTA, // delta_bc*[A*B]_ad
    order_BOTH_OUTER_DELTA    // delta_ad*delta_bc*tr[A*B]
  };

  class Matrix
  {
    std::complex< double > d_data[16];

    public:

      Matrix();
      Matrix(Matrix const &other);
      Matrix(std::complex< double > const &value);
      Matrix(std::complex< double > const * const data);

      std::complex< double > const &operator[](size_t const idx) const;
      std::complex< double >       &operator[](size_t const idx);

      std::complex< double > trace() const;

      Matrix operator+(Matrix const &other) const;
      Matrix operator-(Matrix const &other) const;

      void operator+=(Matrix const &other);
      void operator-=(Matrix const &other);

      bool operator==(Matrix const &rhs);

      template< size_t Index >
      Matrix operator*(Gamma< Index > const &gamma) const;

      template< size_t Index >
      void operator*=(Gamma< Index > const &gamma);
      
      template< size_t Index >
      void right_multiply(Gamma< Index > const &gamma); //applies matrix*=gamma, supplied for consistency

      template< size_t Index >
      void left_multiply(Gamma< Index > const &gamma); //gamma * matrix

      Matrix operator*(double const &factor) const;
      Matrix operator*(std::complex< double > const &factor) const;

      void operator*=(double const &factor);
      void operator*=(std::complex< double > const &factor);

      void operator*=(Matrix const &rhs);
      Matrix operator*(Matrix const &rhs) const;

      void transpose();

      // this is just a simple multiplication of the kind (C)_ij = (A)_ij * (B)_ij (no sum!)
      Matrix elementwise_product(Matrix const &other) const;

      // returns array of 16 Matrices in result
      // void outer_product(Matrix const &other, Matrix* result) const;
      void outer_product(Matrix const &other, std::complex< double > * const result, OuterProductIndexOrder const idxOrd) const;

      // returns the data of a Dirac::Matrix in result. source and sink index are fixed
      void outer_product(Matrix const &other, std::complex< double > * const result, OuterProductIndexOrder const idxOrd,
                         Base::DiracIndex sourceIndex, Base::DiracIndex sinkIndex ) const;

      // needed for threepoints
      void eq_sandwich_operator(Matrix const &first, Base::Operator const op, Matrix const &second);

      std::complex< double > const &operator()(Base::DiracIndex const Dirac_src, Base::DiracIndex const Dirac_snk) const;

      size_t size() const;

      template< size_t Index >
      friend Matrix operator*(Gamma< Index > const &gamma, Matrix const &mat);

      friend std::ostream &operator<<(std::ostream &out, Matrix const &mat);

  };

  template< size_t Index >
  Matrix operator*(Gamma< Index > const &gamma, Matrix const &mat);

  std::ostream &operator<<(std::ostream &out, Matrix const &mat);

  #include "Matrix/Matrix.inlines"
  #include "Matrix/Matrix.constructors.inlines"
  #include "Matrix/Matrix.operators.inlines"
  #include "Matrix/Matrix_operator.gamma.templates"
}
