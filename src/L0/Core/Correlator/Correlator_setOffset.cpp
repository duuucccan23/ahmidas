#include "Correlator.ih"

namespace Core
{
  void Correlator::setOffset(size_t timeslice)
  {
    isolate();
    d_offset = timeslice % T();
  }
}
