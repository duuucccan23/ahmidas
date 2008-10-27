#ifndef GUARD_QCD_TENSOR_H
#define GUARD_QCD_TENSOR_H

#include <complex>

#include <L0/Base/Base.h>

namespace QCD
{
  class Tensor;
  class hcTensor;

  std::complex< double > tr(Tensor const &tensor);
  std::complex< double > tr(hcTensor const &tensor);
  
  class Tensor
  {
    std::complex< double > d_data[144];

    public:
      Tensor();
      Tensor(Tensor const &other);
      Tensor(std::complex< double > *data);
      explicit Tensor(hcTensor const &other);
      Tensor &operator=(Tensor const &other);
      Tensor &operator=(hcTensor const &other);

      std::complex< double > &operator[](size_t const idx);
      std::complex< double > const &operator[](size_t const idx) const;

      std::complex< double > &operator()(size_t const dirSink, size_t const colSink,
                                         size_t const dirSource, size_t const colSource);
      std::complex< double > const &operator()(size_t const dirSink, size_t const colSink,
                                               size_t const dirSource, size_t const colSource) const;

      Tensor &leftMultiply(Tensor const &other);
      Tensor &rightMultiply(Tensor const &other);

      Tensor &leftMultiply(hcTensor const &other);
      Tensor &rightMultiply(hcTensor const &other);

      hcTensor dagger() const;
      std::complex< double > trace() const;

      size_t size() const;
  };

  class hcTensor
  {
    friend class Tensor;

    Tensor const &d_parent;

    hcTensor(Tensor const &parent);

    public:
      hcTensor(hcTensor const &parent);

      std::complex< double > operator[](size_t const idx) const;

      Tensor const &dagger() const;

      size_t size() const;
  };
}

#include "Tensor/Tensor.inlines"
#include "Tensor/hcTensor.inlines"

#endif
