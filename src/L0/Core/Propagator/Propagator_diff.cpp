#include "Propagator.ih"

namespace Core
{
  double Propagator::diff(Propagator const& other) const
  {
    if (d_components == other.d_components)
      return 0.0;
    if (L() != other.L() || T() != other.T())
      return -1.0;

    double diff(0.0);

    Propagator::const_iterator it1 = begin();
    Propagator::const_iterator it2 = other.begin();

    while(it1 != end())
    {
      diff += abs((*it1).diff(*it2));
      ++it2;
      ++it1;
    }
    diff /= (double(d_components->volume())*double(size()));
    return diff;
  }
}
