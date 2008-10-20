#ifndef GUARD_QCD_TENSOR_H
#define GUARD_QCD_TENSOR_H

#include <complex>

#include <L0/Base/Base.h>

namespace QCD
{
  class Tensor
  {
    std::complex< double > d_data[144];

    public:
      Tensor();
      Tensor(Tensor const &other);
      Tensor(std::complex< double > *data);

      std::complex< double > &operator[](size_t const idx);
      std::complex< double > const &operator[](size_t const idx) const;

      std::complex< double > &operator()(size_t const dirSink, size_t const colSink, size_t const dirSource, size_t const colSource);
      std::complex< double > const &operator()(size_t const dirSink, size_t const colSink, size_t const dirSource, size_t const colSource) const;

      Tensor &leftMultiply(Tensor const &other);
      Tensor &rightMultiply(Tensor const &other);
  };
}

#include "Tensor/Tensor.inlines"

#endif
