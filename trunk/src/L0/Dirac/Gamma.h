#pragma once

#include <complex>
#include <iostream>
#include <iomanip>

namespace Dirac
{
  template< size_t Index >
  class Gamma
  {
    static size_t const s_perm[4];
    static std::complex< double > const s_sign[4];

    public:
      size_t const &perm(size_t index) const;
      std::complex< double > const &sign(size_t index) const;
  };
  template< size_t Index >
  std::ostream &operator<<(std::ostream &out, Dirac::Gamma< Index > const &gam);

}

#include "Gamma/Gamma.inlines"
#include "Gamma/Gamma_cout_operator_lshift.template"