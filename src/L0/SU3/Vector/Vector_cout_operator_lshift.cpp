#include "Vector.ih"

std::ostream &SU3::operator<<(std::ostream &out, SU3::GeneralVector< 1 > const &vec)
{
  out << std::scientific << std::setprecision(2) << std::showpos
      << "[ " << vec[0].real() << "  " << vec[0].imag() << " * i   "
              << vec[1].real() << "  " << vec[1].imag() << " * i   "
              << vec[2].real() << "  " << vec[2].imag() << " * i ]"
      	      << std::endl;
  return out;
}
