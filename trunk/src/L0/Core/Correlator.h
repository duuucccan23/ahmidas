#pragma once

#include <algorithm>
#include <numeric>
#include <complex>
#include <string>

#include <L0/Base/Weave.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Tensor.h>


namespace Core
{

  class Correlator
  {

    friend class Weave;

    size_t T;
    size_t L;
    Base::Weave *d_weave;
    size_t *d_references;
    Field < QCD::reducedTensor > *d_data;
    QCD::reducedTensor *d_sumTimeslice;

    // for compatibility with parallel code
    QCD::reducedTensor *d_sumTimeslice_global;

    public:

      Correlator(size_t const L_, size_t const T_, Field < QCD::reducedTensor > *d_data);
      Correlator(Correlator const &other);

      ~Correlator();

      QCD::reducedTensor &operator[](size_t const idx);
      QCD::reducedTensor const &operator[](size_t const idx) const;

      void operator*=(double const factor);
      void operator*=(std::complex< double > const &factor);
      void operator*=(Base::BaryonPropagatorProjector const projector);
      void operator+=(Correlator const &other);

      void sumOverSpatialVolume();
      void sumOverSpatialVolume(size_t const *momentum);

      std::complex <double> getTrSum(size_t const timeslice) const;

      bool isRoot() const;

      size_t getT() const;
      size_t getL() const;
      size_t size() const;


      friend std::ostream &operator<<(std::ostream &out, Correlator const &c);

    private:
      void destroy();
      void isolate();
  };

  std::ostream &operator<<(std::ostream &out, Correlator const &c);


}

#include "Correlator/Correlator.inlines"

#include "Correlator/Correlator_destroy.template"
#include "Correlator/Correlator_isolate.template"
