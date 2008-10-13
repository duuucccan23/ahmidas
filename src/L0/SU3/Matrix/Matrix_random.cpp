#include "Matrix.ih"

SU3::Matrix SU3::Matrix::random()
{
  Matrix result;
  std::generate_n(result.d_data, 18, Base::Random::uniform);
  return result;
}
