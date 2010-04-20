#include "Correlator.ih"

std::ostream &Core::operator<<(std::ostream &out, Core::Correlator const &c)
{
  size_t prec_tmp = std::cout.precision();
  out << "T        real part   imaginary part\n";
  for (size_t t = 0; t < c.T(); t++)
  {
    out << t << "  " << std::scientific << std::showpos << std::setprecision(8);
    out << (c[t]).trace().real() << "  " << (c[t]).trace().imag() << std::endl;
  }
  std::cout.precision(prec_tmp);
  return out;
}
