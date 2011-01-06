#include "Matrix.ih"

// pre
bool SU3::Matrix::equals(SU3::Matrix const &other, double const relativePrecision) const
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
    if (Base::isNaN(abs(d_data[idx])))
      return false;
    if (abs(tmp.d_data[idx]) != ZERO && (abs(tmp.d_data[idx])/abs(d_data[idx])) > relativePrecision)
      return false;
  }
  return true;
}
