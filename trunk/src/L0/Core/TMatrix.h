#ifndef GUARD_CORE_TMATRIX_H
#define GUARD_CORE_TMATRIX_H

#include <algorithm>
#include <L0/QCD/Tensor.h>

namespace Core
{
  template< size_t T >
  class TMatrix
  {
    size_t      *d_references;
    QCD::Tensor *d_transfer;

    public:
      TMatrix();
      TMatrix(QCD::Tensor *transfer);
      TMatrix(TMatrix const &other);
      TMatrix &operator=(TMatrix const &other);
      
      ~TMatrix();

    private:
      void destroy();
      void isolate();
  }
}

#include "TMatrix/TMatrix.inlines"

#include "TMatrix/TMatrix_destroy.template"
#include "TMatrix/TMatrix_isolate.template"
#include "TMatrix/TMatrix_operator_assign.template"

#endif
