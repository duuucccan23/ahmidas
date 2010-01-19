#pragma once

#include <algorithm>
#include <complex>
#include <string>

#include <L0/Base/Weave.h>


namespace Core
{
  class Correlator
  {
    size_t T;
    size_t L;
    Base::Weave            &d_weave;
    size_t                 *d_references;
    std::complex< double > *d_data;

    public:
      Correlator(size_t const L_, size_t const T_);
      Correlator(size_t const L_, size_t const T_, std::complex< double > const &value);
      Correlator(Correlator const &other);

      ~Correlator();

      std::complex< double > &operator[](size_t const idx);
      std::complex< double > const &operator[](size_t const idx) const;

      //void sumOverTimeSlices();

      void save(std::string const&file);

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
