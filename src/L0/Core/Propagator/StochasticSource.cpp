#include "StochasticSource.ih"

namespace Core
{
  Propagator StochasticSource::operator*(Propagator const &propagator) const
  {

    Propagator tmp(propagator);
    assert (T()==propagator.T() && L()==propagator.L());



    Propagator::iterator Itmp(tmp.begin()); //isolate() is called automatically here
    Propagator::const_iterator Is(begin());
    Propagator::const_iterator Ip(begin());

    while(Is != end())
    {
      std::transform(&((*Is)[0]), &((*Is)[0]) + 144, &((*Ip)[0]), &((*Itmp)[0]),
                     std::multiplies< std::complex< double > >());
      ++Itmp;
      ++Is;
      ++Ip;
    }
    assert(Itmp == tmp.end());
    assert(Ip == propagator.end());
    return tmp;
  }

}
