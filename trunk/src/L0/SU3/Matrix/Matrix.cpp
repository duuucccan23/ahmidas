#include "Matrix.ih"

namespace
{
  double const identity[] =
  {
    1.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    1.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    1.0, 0.0
  };

  double const zero[] =
  {
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0,
    0.0, 0.0
  };
}

SU3::Matrix const SU3::Matrix::s_identity(::identity);
SU3::Matrix const SU3::Matrix::s_zero(::zero);

