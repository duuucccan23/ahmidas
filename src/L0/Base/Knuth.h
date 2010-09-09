#pragma once

#include <cstdlib>
#include <stdint.h>

namespace Base
{
  class Knuth
  {
    static Knuth *s_instance;

    size_t d_next;
    size_t d_lag;

    uint64_t const d_max;
    uint64_t const d_seed;

    double d_state[56];

    public:
      static Knuth &instance(uint64_t seed = 0);
      double operator()();

      ~Knuth();

    private:
      Knuth(uint64_t const seed);
      void initialize(uint64_t seed = 0);
  };
}

#include "Knuth/Knuth.inlines"
