#include "Matrix.ih"

namespace
{
  double randDouble()
  {
    return (1 - (2.0 * rand()) / RAND_MAX);
  }
}

SU3::Matrix SU3::Matrix::random(int seed)
{
  if (seed)
    srand(seed);
  else
    srand(time(0));
  double data[18];
  std::generate(data, data + 18, randDouble);
  return Matrix(data);
}
