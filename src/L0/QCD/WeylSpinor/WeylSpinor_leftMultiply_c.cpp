#include "WeylSpinor.ih"

namespace QCD
{ 
  template< >
  void QCD::WeylSpinor::leftMultiply(Dirac::Pauli< 2 >, std::complex< double > const &factor)
  {
    SU3::Vector temp(d_data[0]);
    std::transform(reinterpret_cast< std::complex< double >* >(d_data + 1), 
                  reinterpret_cast< std::complex< double >* >(d_data + 1) + 3,
                  reinterpret_cast< std::complex< double >* >(d_data),
                  std::bind1st(std::multiplies< std::complex< double > >(), factor * std::complex< double >(0.0, -1.0)));
    std::transform(reinterpret_cast< std::complex< double >* >(&temp), 
                  reinterpret_cast< std::complex< double >* >(&temp) + 3,
                  reinterpret_cast< std::complex< double >* >(d_data + 1),
                  std::bind1st(std::multiplies< std::complex< double > >(), factor * std::complex< double >(0.0, 1.0)));
  }
}
