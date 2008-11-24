#include "Vector.ih"

void SU3::Vector::leftMultiply(SU3::hcMatrix const &mat)
{
  SU3::Matrix const &parent = mat.dagger();
  std::complex< double > temp[] = {d_data[0], d_data[1], d_data[2]};
  d_data[0] = std::conj(parent(0, 0)) * temp[0] + std::conj(parent(1, 0)) * temp[1] + std::conj(parent(2, 0)) * temp[2];
  d_data[1] = std::conj(parent(0, 1)) * temp[0] + std::conj(parent(1, 1)) * temp[1] + std::conj(parent(2, 1)) * temp[2];
  d_data[2] = std::conj(parent(0, 2)) * temp[0] + std::conj(parent(1, 2)) * temp[1] + std::conj(parent(2, 2)) * temp[2];
}
