#pragma once

#include <L0/Base/Weave.h>
#include <L0/Dirac/Matrix.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Tensor.h>

namespace Core
{

  class Correlator
  {

    friend class Weave;

    Base::Weave *d_weave;
    size_t *d_references;
    Field < Dirac::Matrix > *d_data;
    Dirac::Matrix *d_sumTimeslice;

    // for compatibility with parallel code
    Dirac::Matrix *d_sumTimeslice_global;

    public:

      Correlator(size_t const L, size_t const T, Field < Dirac::Matrix > *d_data);
      Correlator(Correlator const &other);

      ~Correlator();

      Dirac::Matrix &operator[](size_t const idx);
      Dirac::Matrix const &operator[](size_t const idx) const;

      void operator*=(double const factor);
      void operator*=(std::complex< double > const &factor);
      void operator*=(Base::BaryonPropagatorProjector const projector);
      void operator+=(Correlator const &other);

      void sumOverSpatialVolume();
      void sumOverSpatialVolume(size_t const *momentum);

      std::complex <double> getTrSum(size_t const timeslice) const;

      bool isRoot() const;

      size_t T() const;
      size_t L() const;
      size_t size() const;

      friend std::ostream &operator<<(std::ostream &out, Correlator const &c);

    private:
      void destroy();
      void isolate();
  };

  std::ostream &operator<<(std::ostream &out, Correlator const &c);

}

#include "Correlator/Correlator.inlines"
