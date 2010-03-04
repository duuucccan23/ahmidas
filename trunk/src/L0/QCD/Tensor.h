#pragma once

#include <complex>
#include <algorithm>
#include <functional>

#include <L0/Base/Base.h>
#include <L0/Base/Random.h>
#include <L0/Dirac/Gamma.h>
#include <L0/QCD/Spinor.h>


namespace QCD
{
  class Tensor;
  class hcTensor;
  class reducedTensor;

  std::complex< double > tr(Tensor const &tensor);
  std::complex< double > tr(hcTensor const &tensor);
  std::complex< double > tr(reducedTensor const &rTensor);

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

    // Tensor without colour structure
    friend class reducedTensor;

    std::complex< double > d_data[144];

    //those functions are for construction of baryon fields
    void left_multiply_proton();
    void right_multiply_proton();

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

      /*  NOTE for performance reason the (template) member operator
          Tensor Tensor::operator*(Dirac::Gamma< Index > const &gamma) const
          and the non-member operator
          Tensor operator*(Dirac::Gamma< Index > const &gamma, Tensor const &tensor)
          should be replaced by
          Tensor::left_multiply(Dirac::Gamma< Index > const &gamma, Tensor &result) const
          and
          Tensor::right_multiply(Dirac::Gamma< Index > const &gamma, Tensor &result) const
          in the near future.
      */
      template< size_t Index >
      Tensor operator*(Dirac::Gamma< Index > const &gamma) const;

      template< size_t Index >
      Tensor operator*(Dirac::Sigma< Index > const &gamma) const;

      template< size_t Index >
      void operator*=(Dirac::Gamma< Index > const &gamma);

      template< size_t Index >
      void operator*=(Dirac::Sigma< Index > const &gamma);


      Tensor &leftMultiply(Tensor const &other);
      Tensor &rightMultiply(Tensor const &other);

      Tensor &leftMultiply(hcTensor const &other);
      Tensor &rightMultiply(hcTensor const &other);

      hcTensor dagger() const;
      void transposeDirac();
      std::complex< double > trace() const;

      size_t size() const;
      std::complex< double > diff(Tensor const &other) const;
      void setToRandom();

  #include "Tensor/Tensor.iterator"

      iterator begin(Base::ColourIndex const idx, TensorColourStride const stride);

      iterator end(Base::ColourIndex const idx, TensorColourStride const stride);

      iterator begin(Base::DiracIndex const idx, TensorDiracStride const stride);

      iterator end(Base::DiracIndex const idx, TensorDiracStride const stride);

      friend std::ostream &operator<<(std::ostream &out, Tensor const &tensor);

      template< size_t Index >
      friend Tensor operator*(Dirac::Gamma< Index > const &gamma, Tensor const &tensor);
      template< size_t Index >
      friend Tensor operator*(Dirac::Sigma< Index > const &gamma, Tensor const &tensor);
  };

  std::ostream &operator<<(std::ostream &out, Tensor const &tensor);

  template< size_t Index >
  QCD::Tensor operator*(Dirac::Gamma< Index > const &gamma, Tensor const &tensor);
  template< size_t Index >
  QCD::Tensor operator*(Dirac::Sigma< Index > const &gamma, Tensor const &tensor);


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

  /* reduced tensor a.k.a. "Dirac matrix" */
  class reducedTensor
  {
    friend class Tensor;

    std::complex< double > d_data[16];

    reducedTensor(Tensor const &fullTensor, Base::ColourIndex const colour_src, Base::ColourIndex const colour_snk);

    public:

// old idea:
//       // for mesons: colours are summed over with delta_ab
//       template< size_t Index >
//       reducedTensor(Tensor const &a, Tensor const &b, Dirac::Gamma< Index > const &Dirac_structure);
//       // for baryons: colours are summed over with epsilon_abc (Dirac structure is assumed to involve b and c)
//       template< size_t Index >
//       reducedTensor(Tensor const &a, Tensor const &b, Tensor const &c, Dirac::Gamma< Index > const &Dirac_structure);

      reducedTensor();
      reducedTensor(reducedTensor const &other);
      reducedTensor(std::complex< double > const &value);

      // important for contractions
      reducedTensor(Tensor const &a, Tensor const &b);
      reducedTensor(Tensor const &A, Tensor const &B, bool const colourDilutedSource);
      reducedTensor(Tensor const &A, Tensor const &B, Tensor const &C, Base::BaryonInterpolatingField iPol);

      std::complex< double > trace() const;

      reducedTensor operator+(reducedTensor const &other) const;
      reducedTensor operator-(reducedTensor const &other) const;

      void operator+=(reducedTensor const &other);
      void operator-=(reducedTensor const &other);

      reducedTensor &operator=(reducedTensor const &rhs);

      template< size_t Index >
      reducedTensor operator*(Dirac::Gamma< Index > const &gamma) const;

      template< size_t Index >
      void operator*=(Dirac::Gamma< Index > const &gamma);
      template< size_t Index >
      void left_multiply(Dirac::Gamma< Index > const &gamma); //does the same

      reducedTensor operator*(double const &factor) const;
      reducedTensor operator*(std::complex< double > const &factor) const;

      void operator*=(double const &factor);
      void operator*=(std::complex< double > const &factor);

      void operator*=(reducedTensor const &rhs);
      reducedTensor operator*(reducedTensor const &rhs) const;

      std::complex< double > const &operator()(Base::DiracIndex const Dirac_src, Base::DiracIndex const Dirac_snk) const;

      size_t size() const;

      template< size_t Index >
      friend reducedTensor operator*(Dirac::Gamma< Index > const &gamma, reducedTensor const &rTensor);

      friend std::ostream &operator<<(std::ostream &out, reducedTensor const &rTensor);
  };

  template< size_t Index >
  QCD::reducedTensor operator*(Dirac::Gamma< Index > const &gamma, reducedTensor const &rTensor);

  std::ostream &operator<<(std::ostream &out, reducedTensor const &rTensor);

}

#include "Tensor/Tensor.inlines"
#include "Tensor/hcTensor.inlines"
#include "Tensor/Tensor.gamma.inlines"
#include "Tensor/Tensor.sigma.inlines"

#include "Tensor/reducedTensor.inlines"
#include "Tensor/reducedTensor.constructors.inlines"
#include "Tensor/reducedTensor.operators.inlines"
#include "Tensor/reducedTensor.gamma.inlines"

#include "Tensor/Tensor.iterator.inlines"
