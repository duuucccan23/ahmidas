#include "Vector.ih"

void SU3::Vector::leftMultiply(SU3::Matrix const &mat)
{
  std::complex< double > temp[] = {d_data[0], d_data[1], d_data[2]};
  d_data[0] = mat(0, 0) * temp[0] + mat(0, 1) * temp[1] + mat(0, 2) * temp[2];
  d_data[1] = mat(1, 0) * temp[0] + mat(1, 1) * temp[1] + mat(1, 2) * temp[2];
  d_data[0] = mat(2, 0) * temp[0] + mat(2, 1) * temp[1] + mat(2, 2) * temp[2];
}
