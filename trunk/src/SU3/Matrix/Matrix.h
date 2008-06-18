#ifndef GUARD_SU3_MATRIX_H
#define GUARD_SU3_MATRIX_H

#include <complex>

namespace SU3
{
  class Matrix;
  class hcMatrix;
  
  std::complex< double > det(Matrix const &mat);
  std::complex< double > det(hcMatrix const &mat);
  
  std::complex< double > tr(Matrix const &mat);
  std::complex< double > tr(hcMatrix const &mat);
  
  std::complex< double > minorSum(Matrix const &mat);
  std::complex< double > minorSum(hcMatrix const &mat);
  
  class Matrix
  {
    std::complex< double > d_data[9];
    
    // Some useful constant matrices
    static const Matrix s_identity;
    static const Matrix s_zero;   
    
    public:
      Matrix();
      Matrix(double *data);
      Matrix(std::complex< double > *data);      
      Matrix(Matrix const &other);
      explicit Matrix(hcMatrix const &other);

      Matrix &operator=(Matrix const &other);
      
      static Matrix const &identity();
      static Matrix const &zero();
      
      ~Matrix();
      
      std::complex< double > &operator()(int i, int j);
      std::complex< double > const &operator()(int i, int j) const;

      template< typename T >
      Matrix &operator+=(T const &rhand);
      
      template< typename T >
      Matrix &operator-=(T const &rhand);
      
      template< typename T >
      Matrix &operator*=(T const &rhand);

      template< typename T >
      Matrix &operator/=(T const &rhand);

      
      void leftMultiply(Matrix const &other);
      void rightMultiply(Matrix const &other);
      
      void leftMultiply(hcMatrix const &other);
      void rightMultiply(hcMatrix const &other);           
      
      hcMatrix const dagger();
      
      void reunitarize();
      
      friend std::complex< double > det(Matrix const &mat);
      friend std::complex< double > tr(Matrix const &mat);
      friend std::complex< double > minorSum(Matrix const &mat);
      
    private:
      void givens(std::complex< double > &c, std::complex< double > &s, 
                  std::complex< double > const &f, std::complex< double > const &g);
      std::complex< double > sign(std::complex< double > const &x);
      
      std::complex< double > det(std::complex< double > const *data);
      std::complex< double > tr(std::complex< double > const *data);
      std::complex< double > minorSum(std::complex< double > const *data);
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
      friend std::complex< double > minorSum(hcMatrix const &mat);
  };
}

#include "Matrix.inlines"
#include "hcMatrix.inlines"

#endif
