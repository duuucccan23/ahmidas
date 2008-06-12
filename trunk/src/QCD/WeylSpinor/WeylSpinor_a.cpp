#include "WeylSpinor.ih"

namespace QCD
{
  template< >
  void WeylSpinor::leftMultiply(Dirac::Pauli< 1 >, std::complex< double > const &factor)
  {
    // NOTE A lot of tricky casting going on here... Double check this!
    // I attempt to exploit SU3::Vector's internal structure as a three complex doubles...
    // Lots of clutter can be removed here if we accept a less optimized copy -- consider this.
    SU3::Vector temp(d_data[0]);
    std::transform(reinterpret_cast< std::complex< double >* >(d_data + 1), 
                  reinterpret_cast< std::complex< double >* >(d_data + 1) + 3,
                  reinterpret_cast< std::complex< double >* >(d_data),
                  std::bind1st(std::multiplies< std::complex< double > >(), factor));
    std::transform(reinterpret_cast< std::complex< double >* >(&temp), 
                  reinterpret_cast< std::complex< double >* >(&temp) + 3,
                  reinterpret_cast< std::complex< double >* >(d_data + 1),
                  std::bind1st(std::multiplies< std::complex< double > >(), factor));
  }
}
