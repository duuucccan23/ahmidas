#include "Matrix.ih"

SU3::Matrix &SU3::Matrix::operator=(SU3::hcMatrix const &other)
{
  if (&(other.d_parent) != this)
  {
    d_data[0] = std::conj(other.d_parent.d_data[0]);
    d_data[1] = std::conj(other.d_parent.d_data[3]);
    d_data[2] = std::conj(other.d_parent.d_data[6]);
    d_data[3] = std::conj(other.d_parent.d_data[1]);
    d_data[4] = std::conj(other.d_parent.d_data[4]);
    d_data[5] = std::conj(other.d_parent.d_data[7]);
    d_data[6] = std::conj(other.d_parent.d_data[2]);
    d_data[7] = std::conj(other.d_parent.d_data[5]);
    d_data[8] = std::conj(other.d_parent.d_data[8]);
  }
  else
  {
    std::complex< double > temp[3] = {d_data[1], d_data[2], d_data[5]};
    d_data[0] = std::conj(d_data[0]);
    d_data[1] = std::conj(d_data[3]);
    d_data[2] = std::conj(d_data[6]);
    d_data[3] = std::conj(temp[0]);
    d_data[4] = std::conj(d_data[4]);
    d_data[5] = std::conj(d_data[7]);
    d_data[6] = std::conj(temp[1]);
    d_data[7] = std::conj(temp[2]);
    d_data[8] = std::conj(d_data[8]);
  }
  return *this;
}
