#include "Matrix.ih"

namespace
{
  double const ar_identity[] =
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

  double const ar_zero[] =
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

SU3::Matrix const SU3::Matrix::s_identity(ar_identity);
SU3::Matrix const SU3::Matrix::s_zero(ar_zero);

