#include "Matrix.ih"

void SU3::Matrix::leftMultiply(SU3::Matrix const &other)
{
  std::complex< double > temp[9];
  
  temp[0] = other.d_data[0] * d_data[0] + other.d_data[1] * d_data[3] + other.d_data[2] * d_data[6];
  temp[1] = other.d_data[0] * d_data[1] + other.d_data[1] * d_data[4] + other.d_data[2] * d_data[7];
  temp[2] = other.d_data[0] * d_data[2] + other.d_data[1] * d_data[5] + other.d_data[2] * d_data[8];
  temp[3] = other.d_data[3] * d_data[0] + other.d_data[4] * d_data[3] + other.d_data[5] * d_data[6];
  temp[4] = other.d_data[3] * d_data[1] + other.d_data[4] * d_data[4] + other.d_data[5] * d_data[7];
  temp[5] = other.d_data[3] * d_data[2] + other.d_data[4] * d_data[5] + other.d_data[5] * d_data[8];
  temp[6] = other.d_data[6] * d_data[0] + other.d_data[7] * d_data[3] + other.d_data[8] * d_data[6];
  temp[7] = other.d_data[6] * d_data[1] + other.d_data[7] * d_data[4] + other.d_data[8] * d_data[7];
  temp[8] = other.d_data[6] * d_data[2] + other.d_data[7] * d_data[5] + other.d_data[8] * d_data[8];

  std::copy(temp, temp + 9, d_data);
}

