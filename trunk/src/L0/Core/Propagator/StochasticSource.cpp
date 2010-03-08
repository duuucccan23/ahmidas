#include "Propagator.ih"

namespace Core
{
  template<  >
  Propagator StochasticSource< 4 >::operator*(StochasticPropagator< 4 > const &sPropagator) const
  {

    Propagator tmp(sPropagator);
    assert (T()==sPropagator.T() && L()==sPropagator.L());

//     if (!(NComp == 4))
//     {
//       std::cerr << "Propagator StochasticSource< NComp >::operator*(StochasticPropagator< NComp > const &) const\n"
//           << "has not been implemented yet for NComp != 4" << std::endl;
//       exit(1);
//     }

    Propagator::iterator Itmp(tmp.begin()); //isolate() is called automatically here
    Propagator::const_iterator Is(begin());
    Propagator::const_iterator Ip(sPropagator.begin());

    while(Is != end())
    {
      std::transform(&((*Is)[0]), &((*Is)[0]) + 144, &((*Ip)[0]), &((*Itmp)[0]),
                     std::multiplies< std::complex< double > >());
      ++Itmp;
      ++Is;
      ++Ip;
    }
    assert(Itmp == tmp.end());
    assert(Ip == sPropagator.end());
    return tmp;
  }

}
