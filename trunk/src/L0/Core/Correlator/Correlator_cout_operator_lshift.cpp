#include "Correlator.ih"

namespace Core
{

  template< >
  std::ostream &operator<<(std::ostream &out, Correlator< std::complex< double > > const &c)
  {
    if(!c.d_weave->isRoot())
      return out;

    size_t prec_tmp = std::cout.precision();
    out << "T        real part   imaginary part\n";
    for (size_t t = 0; t < c.T(); t++)
    {
      out.width(3);
      out << t << "  " << std::scientific << std::showpos << std::setprecision(8);
      out << (c[t]).real() << "  " << (c[t]).imag() << std::endl;
    }
    std::cout.precision(prec_tmp);
    return out;
  }

  template< >
  std::ostream &operator<<(std::ostream &out, Correlator< Dirac::Matrix > const &c)
  {
    if(!c.d_weave->isRoot())
      return out;

    size_t prec_tmp = std::cout.precision();
    out << "T        real part   imaginary part\n";
    for (size_t t = 0; t < c.T(); t++)
    {
      out.width(3);
      out << t << "  " << std::scientific << std::showpos << std::setprecision(8);
      out << (c[t]).trace().real() << "  " << (c[t]).trace().imag() << std::endl;
    }
    std::cout.precision(prec_tmp);
    return out;
  }
}
