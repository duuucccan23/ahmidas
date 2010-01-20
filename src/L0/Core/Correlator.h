#pragma once

#include <algorithm>
#include <complex>
#include <string>

#include <L0/Base/Weave.h>
#include <L0/Core/Field.h>
#include <L0/QCD/Tensor.h>

namespace Core
{

  class Correlator
  {
    size_t T;
    size_t L;
    Base::Weave  &d_weave;
    size_t *d_references;
    QCD::reducedTensor *d_data;

    public:
      Correlator(size_t const L_, size_t const T_);
      Correlator(size_t const L_, size_t const T_, std::complex< double > const &value);
      Correlator(Correlator const &other);

      ~Correlator();

      QCD::reducedTensor &operator[](size_t const idx);
      QCD::reducedTensor const &operator[](size_t const idx) const;

      void sumOverTimeSlices();

      void sumOverTimeSlices(Core::Field< QCD::reducedTensor > **timeslices);
      void sumOverTimeSlices(Core::Field< QCD::reducedTensor > **timeslices, double const *momentum);

      void save(std::string const&file);

      QCD::reducedTensor sum(Core::Field< QCD::reducedTensor > const *field);
      QCD::reducedTensor sum(Core::Field< QCD::reducedTensor > const *field, double * const momentum);


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
