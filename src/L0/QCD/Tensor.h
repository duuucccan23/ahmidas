#pragma once

#include <complex>
#include <algorithm>
#include <functional>

#include <L0/Base/Base.h>
#include <L0/Base/Random.h>
#include <L0/Base/Z2.h>
#include <L0/Dirac/Gamma.h>
#include <L0/Dirac/Matrix.h>
#include <L0/QCD/Spinor.h>

namespace QCD
{
  class Tensor;
  class hcTensor;

  std::complex< double > tr(Tensor const &tensor);
  std::complex< double > tr(hcTensor const &tensor);
  std::complex< double > tr(Dirac::Matrix const &rTensor);

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
    friend class Dirac::Matrix;

    std::complex< double > d_data[144];

     // Some useful constant matrices
    static const Tensor s_identity;

    public:
      Tensor();
      Tensor(Tensor const &other);
      Tensor(Spinor *data[12]);
      Tensor(std::complex< double > *data);
      Tensor(Dirac::Matrix const * const data[9]);
      explicit Tensor(hcTensor const &other);

      static Tensor const &identity();

      Tensor &operator=(Tensor const &other);
      Tensor &operator=(hcTensor const &other);

      Spinor &operator[](short const component);
      Spinor const &operator[](short const component) const;

      std::complex< double > &operator()(size_t const idx);
      std::complex< double > const &operator()(size_t const idx) const;

      std::complex< double > &operator()(size_t const dirSink, size_t const colSink,
                                         size_t const dirSource, size_t const colSource);
      std::complex< double > const &operator()(size_t const dirSink, size_t const colSink,
                                               size_t const dirSource, size_t const colSource) const;

      // get entries for fixed color at source and sink
      void getDiracMatrix(Dirac::Matrix &dMatrix,
                          Base::ColourIndex const colour_src,
                          Base::ColourIndex const colour_snk) const;

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
      void rightMultiply(Dirac::Gamma< Index > const &gamma);

      // this is not a multiplication in the sense of the others
      // assumption: xi is diagonal (i.e. spin & color diluted)
      // multiplies diagonal elements of xi to corresponding source entries of *this
      void leftMultiplySpinColorDilutedConj(Tensor const& xi);


      template< size_t Index >
      Tensor operator*(Dirac::Gamma< Index > const &gamma) const;

      Tensor operator*(Dirac::Identity const &gamma) const;

      template< size_t Index >
      Tensor operator*(Dirac::Sigma< Index > const &gamma) const;

      template< size_t Index >
      void operator*=(Dirac::Gamma< Index > const &gamma);

      void operator*=(Dirac::Identity const &gamma);

      template< size_t Index >
      void operator*=(Dirac::Sigma< Index > const &gamma);

      void operator*=(std::complex< double > const &factor);

      void operator+=(Tensor const &other);
      void operator-=(Tensor const &other);

      Tensor &leftMultiply(Tensor const &other);
      Tensor &rightMultiply(Tensor const &other);

      Tensor &leftMultiply(hcTensor const &other);
      Tensor &rightMultiply(hcTensor const &other);

      Tensor &leftMultiply(SU3::Matrix const &mat);
      Tensor &rightMultiply(SU3::Matrix const &mat);

      Tensor &leftMultiply(SU3::hcMatrix const &mat);
      Tensor &rightMultiply(SU3::hcMatrix const &mat);

      //those functions are for construction of baryon fields
      void left_multiply_proton();
      void right_multiply_proton();

      // NOTE: there should be a discussion how this is to be named
      // this function contracts color indices according to epsilon tensor
      // and Dirac indices with delta
      void make_sequential(Tensor const &A, Tensor const &B);


      hcTensor dagger() const;
      const Tensor &spinDilutedDagger() const;
      void spinDilutedConjugate();
      void undilutedConjugate();
      void spinDilutedTransposeFull();

      void conjugate();
      void transposeDirac();
      void transposeFull();
      std::complex< double > trace() const;

      size_t size() const;
      std::complex< double > diff(Tensor const &other) const;
      double norm() const;
      double normq() const;
      void setToRandom();
      // this produces tensor filled with random Z(4) elements
      void setToRandom(Base::SourcePolarization const, Base::SourceColorState const,
                       Base::SourceStochasticTypeFlag const type, uint64_t seed = 0);

  #include "Tensor/Tensor.iterator"

      iterator begin(Base::ColourIndex const idx, TensorColourStride const stride);

      iterator end(Base::ColourIndex const idx, TensorColourStride const stride);

      iterator begin(Base::DiracIndex const idx, TensorDiracStride const stride);

      iterator end(Base::DiracIndex const idx, TensorDiracStride const stride);

      friend std::ostream &operator<<(std::ostream &out, Tensor const &tensor);

      template< size_t Index >
      friend Tensor operator*(Dirac::Gamma< Index > const &gamma, Tensor const &tensor);
      friend Tensor operator*(Dirac::Identity const &gamma, Tensor const &tensor);
      template< size_t Index >
      friend Tensor operator*(Dirac::Sigma< Index > const &gamma, Tensor const &tensor);

      friend void getDiracMatrix(Dirac::Matrix &dMatrix, Tensor const &A, Tensor const &B);
      friend void getSpinDilutedDiracMatrix(Dirac::Matrix &dMatrix, Tensor const &A, Tensor const &B);
      friend void getUndilutedDiracMatrix(Dirac::Matrix &dMatrix, Tensor const &A, Tensor const &B);
      friend void getDiracMatrix(Dirac::Matrix &dMatrix, Tensor const &A, Tensor const &B, Tensor const &C, Base::BaryonInterpolatingField const iPol);
      friend void getDiracMatrix_alternative(Dirac::Matrix &dMatrix, Tensor const &A, Tensor const &B, Tensor const &C, Base::BaryonInterpolatingField const iPol);

      // this is a very unelegant solution needed for threepoints
      friend void multiplyOuterDiracIndices(Tensor const &left, Tensor const * const right, Tensor * const result);

      friend void make_sequential_d(Tensor &result, Tensor const &U1, Tensor const &U2, Base::BaryonPropagatorProjector const projector);
      friend void make_sequential_u(Tensor &result, Tensor const &D,  Tensor const &U,  Base::BaryonPropagatorProjector const projector);
      friend void make_sequential_d(Tensor result[16], Tensor const &U1, Tensor const &U2);
      friend void make_sequential_u(Tensor result[16], Tensor const &D, Tensor const &U);
  };

  std::ostream &operator<<(std::ostream &out, Tensor const &tensor);

  template< size_t Index >
  QCD::Tensor operator*(Dirac::Gamma< Index > const &gamma, Tensor const &tensor);
  QCD::Tensor operator*(Dirac::Identity const &gamma, Tensor const &tensor);
  template< size_t Index >
  QCD::Tensor operator*(Dirac::Sigma< Index > const &gamma, Tensor const &tensor);

  void getDiracMatrix(Dirac::Matrix &dMatrix, Tensor const &A, Tensor const &B);
  void getSpinDilutedDiracMatrix(Dirac::Matrix &dMatrix, Tensor const &A, Tensor const &B);
  void getUndilutedDiracMatrix(Dirac::Matrix &dMatrix, Tensor const &A, Tensor const &B);
  void getDiracMatrix(Dirac::Matrix &dMatrix, Tensor const &A, Tensor const &B, Tensor const &C, Base::BaryonInterpolatingField const iPol);
  // alternative version, gives the same result as function above if Tensors A and C are equal
  // note: one should worry about the correct convention here. Still, both functions give the same result for a proton tow point,
  // provided that A and C are equal (i.e. we use the same propagators for both u)
  void getDiracMatrix_alternative(Dirac::Matrix &dMatrix, Tensor const &A, Tensor const &B, Tensor const &C, Base::BaryonInterpolatingField const iPol);

  void multiplyOuterDiracIndices(Tensor const &left, Tensor const * const right, Tensor * const result);

  void make_sequential_d(Tensor &result, Tensor const &U1, Tensor const &U2, Base::BaryonPropagatorProjector const projector);
  void make_sequential_u(Tensor &result, Tensor const &D,  Tensor const &U1, Base::BaryonPropagatorProjector const projector);
  void make_sequential_d(Tensor result[16], Tensor const &U1, Tensor const &U2);
  void make_sequential_u(Tensor result[16], Tensor const &D, Tensor const &U);

  class hcTensor
  {
    friend class Tensor;

    Tensor const &d_parent;

    hcTensor(Tensor const &parent);

    public:
      hcTensor(hcTensor const &parent);

      Spinor operator[](size_t const idx) const;

      std::complex< double > operator()(size_t const idx) const;

      Tensor const &dagger() const;

      size_t size() const;
  };

}

#include "Tensor/Tensor.inlines"
#include "Tensor/Tensor.operators.inlines"
#include "Tensor/hcTensor.inlines"
#include "Tensor/Tensor.gamma.inlines"
#include "Tensor/Tensor.sigma.inlines"

#include "Tensor/Tensor.iterator.inlines"
