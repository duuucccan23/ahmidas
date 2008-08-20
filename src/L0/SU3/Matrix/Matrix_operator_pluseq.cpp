#include "Matrix.ih"

namespace SU3
{
  template< >
  Matrix &Matrix::operator+=(hcMatrix const &rhand)
  {
    d_data[0] += std::conj(rhand.d_parent.d_data[0]);
    d_data[1] += std::conj(rhand.d_parent.d_data[3]);
    d_data[2] += std::conj(rhand.d_parent.d_data[6]);
    d_data[3] += std::conj(rhand.d_parent.d_data[1]);
    d_data[4] += std::conj(rhand.d_parent.d_data[4]);
    d_data[5] += std::conj(rhand.d_parent.d_data[7]);
    d_data[6] += std::conj(rhand.d_parent.d_data[2]);
    d_data[7] += std::conj(rhand.d_parent.d_data[4]);
    d_data[8] += std::conj(rhand.d_parent.d_data[8]);
    return *this;
  }
}
