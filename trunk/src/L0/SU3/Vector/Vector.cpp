#include "Vector.ih"

namespace
{
  double const zero[] =
  {
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0
  };
}

SU3::Tensor< 1 > const SU3::Tensor< 1 >::s_zero(::zero);
