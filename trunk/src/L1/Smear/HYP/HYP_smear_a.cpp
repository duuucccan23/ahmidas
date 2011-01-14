#include "HYP.ih"

namespace Smear
{
  void HYP::smear(Core::Field< QCD::Gauge> &field, size_t iterations) const
  {
    for (size_t ctr = 0; ctr < iterations; ++ctr)
      smear(field);
  }
}
