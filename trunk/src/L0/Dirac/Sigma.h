#pragma once

// This is a spectre-class, that can only be sensibly used in conjunction with
// template functions deriving their implementation from the specific instantiation
// of Sigma< XXX >.

namespace Dirac
{
  template< size_t Index >
  class Sigma
  {};
}
