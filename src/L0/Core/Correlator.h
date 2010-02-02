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

    public:
      Correlator();
      Correlator(size_t const L_, size_t const T_);
      Correlator(size_t const L_, size_t const T_, Field < QCD::reducedTensor > *d_data);
      Correlator(Correlator const &other);

      ~Correlator();

      QCD::reducedTensor &operator[](size_t const idx);
      QCD::reducedTensor const &operator[](size_t const idx) const;

      void sumOverTimeSlices();
      void sumOverTimeSlices(size_t const *momentum);

      void save(std::string const&file);

      std::complex <double> getTrSum(size_t const timeslice);


      size_t getT() const;
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

// #include "Correlator/Correlator_sumOverTimeSlices.template"
