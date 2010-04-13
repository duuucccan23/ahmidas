#include "Correlator.ih"

std::ostream &Core::operator<<(std::ostream &out, Core::Correlator const &c)
{
  for (size_t t = 0; t < c.T(); t++)
  {
    // this is the way formatted output works in C++
    out.width(3);
    out << t << "  ";
    out.width(20);
    // since the Correlator stores complex numbers, we can access the real and imaginary parts using
    // the complex class member functions real() and imag()
    out << std::fixed << std::setprecision(10) << std::showpos << (c[t]).trace().real() << "  ";
    out.width(20);
    out << std::fixed << std::setprecision(10) << std::showpos << (c[t]).trace().imag() << std::endl;
  }
  return out;
}
