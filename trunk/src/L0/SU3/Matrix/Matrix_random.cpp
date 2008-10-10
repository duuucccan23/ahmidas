#include "Matrix.ih"

SU3::Matrix SU3::Matrix::random()
{
  double data[18];
  std::generate_n(data, 18, Base::Random::uniform);
  return Matrix(data);
}
