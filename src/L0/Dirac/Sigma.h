#ifndef GUARD_DIRAC_SPINOR_H
#define GUARD_DIRAC_SPINOR_H

// This is a spectre-class, that can only be sensibly used in conjunction with
// template functions deriving their implementation from the specific instantiation
// of Sigma< XXX >.

namespace Dirac
{
  template< size_t Index >
  class Sigma
  {};
}

#endif
