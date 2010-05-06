#include "Spinor.ih"

bool QCD::Spinor::equals(QCD::Spinor const &other, double const relPrecision) const
{
  std::complex< double > const * ptr_this  = reinterpret_cast< std::complex< double > const * >(d_data);
  std::complex< double > const * ptr_other = reinterpret_cast< std::complex< double > const * >(other.d_data);

  // this is the best that can happen: matrices are absolutely equal
  if (std::equal(ptr_this, ptr_this+12, ptr_other))
    return true;
  // check whether matrices are equal up to precision
  QCD::Spinor tmp(*this);
  std::complex< double > const ZERO(0, 0);
  tmp -= other;
  std::complex< double > const * ptr_tmp  = reinterpret_cast< std::complex< double > * >(tmp.d_data);

  for (size_t idx=0; idx<12; idx++)
  {
    if (abs(ptr_tmp[idx]) != ZERO && (abs(ptr_tmp[idx])/abs(ptr_this[idx])) > relPrecision)
      return false;
  }
  return true;
}
