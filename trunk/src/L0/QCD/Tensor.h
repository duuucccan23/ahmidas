#ifndef GUARD_QCD_TENSOR_H
#define GUARD_QCD_TENSOR_H

#include <complex>
#include <L0/Core/Core.h>

namespace QCD
{
  class Tensor
  {
    size_t                 *d_references;
    std::complex< double > *d_data;

    public:
      Tensor();
      Tensor(Tensor const &other);
      Tensor(std::complex< double > *data);
  };
}

#include "Tensor/Tensor.inlines"

#endif
