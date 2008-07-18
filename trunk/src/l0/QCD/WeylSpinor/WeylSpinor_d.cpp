#include "WeylSpinor.ih"

void QCD::WeylSpinor::swap(std::complex< double > const &myFac, WeylSpinor other, std::complex< double > const &hisFac)
{
  SU3::Vector temp[] = {d_data[0], d_data[1]};    
  std::transform(reinterpret_cast< std::complex< double >* >(other.d_data), 
                 reinterpret_cast< std::complex< double >* >(other.d_data) + 6,
                 reinterpret_cast< std::complex< double >* >(d_data),
                 std::bind1st(std::multiplies< std::complex< double > >(), hisFac));
  std::transform(reinterpret_cast< std::complex< double >* >(&temp), 
                 reinterpret_cast< std::complex< double >* >(&temp) + 6,
                 reinterpret_cast< std::complex< double >* >(other.d_data),
                 std::bind1st(std::multiplies< std::complex< double > >(), myFac));
}
