#ifndef GUARD_CORE_TMATRIX_H
#define GUARD_CORE_TMATRIX_H

#include <algorithm>
#include <L0/QCD/Tensor.h>

namespace Core
{
  template< size_t L, size_t T >
  class hcTMatrix;
  
  template< size_t L, size_t T >
  class TMatrix
  {
    size_t      *d_references;
    QCD::Tensor *d_transfer;

    public:
      TMatrix();
      TMatrix(QCD::Tensor *transfer);
      TMatrix(TMatrix< L, T > const &other);
      explicit TMatrix(hcTMatrix< L, T > const &other);
      TMatrix &operator=(TMatrix< L, T > const &other);
      TMatrix &operator=(hcTMatrix< L, T > const &other);

      ~TMatrix();

      QCD::Tensor &operator[](size_t const idx);
      QCD::Tensor const &operator[](size_t const idx) const;
      
      hcTMatrix< L, T > dagger() const;
      
    private:
      void destroy();
      void isolate();
  };

  template< size_t L, size_t T >
  class hcTMatrix
  {
    friend class TMatrix< L, T >;
    
    TMatrix< L, T > const &d_parent;

    hcTMatrix(TMatrix< L, T > const &other);

    public:
      hcTMatrix(hcTMatrix< L, T > const &other);

      TMatrix< L, T > const &dagger() const;

      QCD::hcTensor operator[](size_t const idx) const;
  };
}

#include "TMatrix/TMatrix.inlines"
#include "TMatrix/hcTMatrix.inlines"

#include "TMatrix/TMatrix_destroy.template"
#include "TMatrix/TMatrix_isolate.template"
#include "TMatrix/TMatrix_operator_assign_a.template"
#include "TMatrix/TMatrix_operator_assign_b.template"

#endif
