#include "Matrix.ih"
#include <cfloat>

bool SU3::Matrix::operator==(SU3::Matrix const &other) const
{
  // this is the best that can happen: matrices are absolutely equal
  if (std::equal(d_data, d_data+9, other.d_data))
    return true;
  // check whether matrices are equal up to precision
  Matrix tmp(*this);
  std::complex< double > const ZERO(0, 0);
  tmp -= other;
  for (size_t idx=0; idx<9; idx++)
  {
    if (abs(tmp.d_data[idx]) != ZERO && (abs(tmp.d_data[idx])/abs(d_data[idx])) > DBL_EPSILON)
      return false;
  }
  return true;
}
