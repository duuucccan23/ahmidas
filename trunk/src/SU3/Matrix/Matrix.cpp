#include "Matrix.ih"

namespace
{
  double ar_identity[] =
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
  
  double ar_zero[] =
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

const SU3::Matrix SU3::Matrix::s_identity(ar_identity);
const SU3::Matrix SU3::Matrix::s_zero(ar_zero);

