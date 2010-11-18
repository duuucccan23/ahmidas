#include "Correlator.ih"

namespace Core
{

  // this has only effect on d_sumTimeslice (maybe this has to be changed in the future;
  template< >
  void Correlator< Dirac::Matrix >::rotateToPhysicalBasis()
  {
    assert(d_sumTimeslice_global != NULL);
    Dirac::Gamma< 5 > gamma_5;
    std::complex< double > const I(0,1);
    isolate();
    for(size_t t = 0; t < T(); t++)
    {
      Dirac::Matrix second = gamma_5 * d_sumTimeslice_global[t];
      second *= I;
      d_sumTimeslice_global[t] += second;
      second = d_sumTimeslice_global[t] * gamma_5;
      second *= I;
      d_sumTimeslice_global[t] += second;
      d_sumTimeslice_global[t] *= 0.5;
    }
  }

}