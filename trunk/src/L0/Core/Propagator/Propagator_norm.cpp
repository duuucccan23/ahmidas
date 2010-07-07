#include "Propagator.ih"

namespace Core
{
  double Propagator::norm() const
  {
    Base::Weave weave(L(), T());

    double local_result(0);

    for (Propagator::const_iterator idx(begin()); idx!=end(); ++idx)
    {
      local_result += (*idx).norm();
    }

    double result;
    weave.allReduce(&local_result, &result);

    // a QCD::Tensor is made of twelve (QCD::Spinor)s
    result /= 12.0;

    // result /= double(weave.globalVolume());

    return result;
  }
}
