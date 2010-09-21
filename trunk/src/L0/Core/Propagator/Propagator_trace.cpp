#include "Propagator.ih"

namespace Core
{
  Field < std::complex < double > > Propagator::trace() const
  {
    Base::Weave weave(L(), T());
    Field < std::complex < double > > traces(L(), T());

    Field < std::complex < double > >::iterator I_tr(traces.begin());

    for (Propagator::const_iterator idx(begin()); idx!=end(); ++idx)
    {
      (*I_tr) = (*idx).trace();
    }
    return traces;
  }
}
