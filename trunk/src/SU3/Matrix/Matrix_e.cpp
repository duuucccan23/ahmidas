#include "Matrix.ih"

void SU3::Matrix::rightMultiply(SU3::hcMatrix const &other)
{
  std::complex< double > const *orig(other.d_parent.d_data);
  std::complex< double > temp[9];
  
  temp[0] = d_data[0] * std::conj(orig[0]) + d_data[1] * std::conj(orig[1]) + d_data[2] * std::conj(orig[2]);
  temp[1] = d_data[0] * std::conj(orig[3]) + d_data[1] * std::conj(orig[4]) + d_data[2] * std::conj(orig[5]);
  temp[2] = d_data[0] * std::conj(orig[6]) + d_data[1] * std::conj(orig[7]) + d_data[2] * std::conj(orig[8]);
  temp[3] = d_data[3] * std::conj(orig[0]) + d_data[4] * std::conj(orig[1]) + d_data[5] * std::conj(orig[6]);
  temp[4] = d_data[3] * std::conj(orig[3]) + d_data[4] * std::conj(orig[4]) + d_data[5] * std::conj(orig[5]);
  temp[5] = d_data[3] * std::conj(orig[6]) + d_data[4] * std::conj(orig[7]) + d_data[5] * std::conj(orig[8]);
  temp[6] = d_data[6] * std::conj(orig[0]) + d_data[7] * std::conj(orig[1]) + d_data[8] * std::conj(orig[2]);
  temp[7] = d_data[6] * std::conj(orig[3]) + d_data[7] * std::conj(orig[4]) + d_data[8] * std::conj(orig[5]);
  temp[8] = d_data[6] * std::conj(orig[6]) + d_data[7] * std::conj(orig[7]) + d_data[8] * std::conj(orig[8]);

  std::copy(temp, temp + 9, d_data);
}
