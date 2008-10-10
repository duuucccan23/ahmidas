#ifndef GUARD_SU3_MATRIX_H
#define GUARD_SU3_MATRIX_H

#include <complex>
#include <iostream>

namespace SU3
{
  class Matrix;
  class hcMatrix;

  std::complex< double > det(Matrix const &mat);
  std::complex< double > det(hcMatrix const &mat);

  std::complex< double > tr(Matrix const &mat);
  std::complex< double > tr(hcMatrix const &mat);

  std::ostream &operator<<(std::ostream &out, Matrix const &mat);

  class Matrix
  {
    std::complex< double > d_data[9];

    // Some useful constant matrices
    static const Matrix s_identity;
    static const Matrix s_zero;

    public:
      Matrix();
      Matrix(double const *data);
      Matrix(std::complex< double > const *data);
      Matrix(Matrix const &other);
      explicit Matrix(hcMatrix const &other);

      Matrix &operator=(Matrix const &other);

      static Matrix const &identity();
      static Matrix const &zero();
      static Matrix random();

      ~Matrix();

      std::complex< double > &operator()(int i, int j);
      std::complex< double > const &operator()(int i, int j) const;

#include "Matrix/Matrix.operators"

      hcMatrix const dagger() const;

      void reunitarize();

      friend std::complex< double > det(Matrix const &mat);
      friend std::complex< double > tr(Matrix const &mat);
      friend std::ostream &operator<<(std::ostream &out, Matrix const &mat);

    private:
      void givens(std::complex< double > &c, std::complex< double > &s,
                  std::complex< double > const &f, std::complex< double > const &g) const;
      std::complex< double > sign(std::complex< double > const &x) const;

      std::complex< double > det(std::complex< double > const *data) const;
      std::complex< double > tr(std::complex< double > const *data) const;
  };

  class hcMatrix
  {
    friend class Matrix;

    Matrix const &d_parent;

    hcMatrix(Matrix const &parent);

    public:
      Matrix const &dagger() const;

      friend std::complex< double > det(hcMatrix const &mat);
      friend std::complex< double > tr(hcMatrix const &mat);
  };
}

#include "Matrix/Matrix.inlines"
#include "Matrix/hcMatrix.inlines"

#endif
