#include "Gauge.ih"

#include <L0/Base/Base.h>

namespace QCD
{
  std::ostream &operator<<(std::ostream &out, Gauge const &gauge)
  {
    out << " mu = ";
    out.width(3);
    out << std::noshowpos << Base::idx_T << " (Base::idx_T)" << std::endl;
    out << gauge[Base::idx_T];
    out << " mu = ";
    out.width(3);
    out << std::noshowpos << Base::idx_X << " (Base::idx_X)" << std::endl;
    out << gauge[Base::idx_X];
    out << " mu = ";
    out.width(3);
    out << std::noshowpos << Base::idx_Y << " (Base::idx_Y)" << std::endl;
    out << gauge[Base::idx_Y];
    out << " mu = ";
    out.width(3);
    out << std::noshowpos <<  Base::idx_Z << " (Base::idx_Z)" << std::endl;
    out << gauge[Base::idx_Z];
    return out;
  }
}
