#include "Propagator.ih"

namespace Core
{
 Propagator &Propagator::revert()
  {
    Dirac::Gamma<5> gamma5 = Dirac::Gamma<5>();
    isolate();
    Propagator::iterator it = begin();
    while(it != end())
    {
      QCD::Tensor tmp((*it).dagger());
      tmp *= gamma5;
      (*it) = gamma5*tmp;
      ++it;
    }
    return *this;
  }
}
