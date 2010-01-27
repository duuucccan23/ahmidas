#pragma once

#include <complex>

namespace Dirac
{
  template< size_t Index >
  class Gamma
  {
    static size_t const s_perm[4];
    static std::complex< double > const s_sign[4];

    public:
      size_t const &perm(size_t index);
      std::complex< double > const &sign(size_t index);
  };

  typedef Gamma< -1 > Unity;
  typedef Gamma< 4 > Gamma0;
  typedef Gamma< 5 > Gamma5;

}

#include "Gamma/Gamma.inlines"
