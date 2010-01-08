#pragma once

#include <complex>

#include <L0/Base/Base.h>
#include "Spinor.h"

namespace QCD
{
  class Tensor;
  class hcTensor;

  std::complex< double > tr(Tensor const &tensor);
  std::complex< double > tr(hcTensor const &tensor);
  
  enum TensorColourStride
  {
    ColourStrideSink   =  1,
    ColourStrideSource = 12
  };

  enum TensorDiracStride
  {
    DiracStrideSink    =  3,
    DiracStrideSource  = 36
  };

  class Tensor
  {
    friend class Spinor;
    std::complex< double > d_data[144];

    public:
      Tensor();
      Tensor(Tensor const &other);
      Tensor(Spinor *data[12]);
      Tensor(std::complex< double > *data);
      explicit Tensor(hcTensor const &other);
      Tensor &operator=(Tensor const &other);
      Tensor &operator=(hcTensor const &other);

      Spinor &operator()(size_t const idx);
      Spinor const &operator()(size_t const idx) const;
      
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
          

  #include "Tensor/Tensor.iterator"

      iterator begin(Base::ColourIndex const idx, TensorColourStride const stride);

      iterator end(Base::ColourIndex const idx, TensorColourStride const stride);

      iterator begin(Base::DiracIndex const idx, TensorDiracStride const stride);

      iterator end(Base::DiracIndex const idx, TensorDiracStride const stride);

  };

  class hcTensor
  {
    friend class Tensor;

    Tensor const &d_parent;

    hcTensor(Tensor const &parent);

    public:
      hcTensor(hcTensor const &parent);

      Spinor operator()(size_t const idx) const;
      
      std::complex< double > operator[](size_t const idx) const;

      Tensor const &dagger() const;

      size_t size() const;
  };
}

#include "Tensor/Tensor.inlines"
#include "Tensor/hcTensor.inlines"

#include "Tensor/Tensor.iterator.inlines"
