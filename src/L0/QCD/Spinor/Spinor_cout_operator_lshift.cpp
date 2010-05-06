#include "Spinor.ih"

std::ostream &QCD::operator<<(std::ostream &out, QCD::Spinor const &spinor)
{
  out << spinor.d_data[0] << spinor.d_data[1] << spinor.d_data[2] << spinor.d_data[3] << std::endl;
  return out;
}

