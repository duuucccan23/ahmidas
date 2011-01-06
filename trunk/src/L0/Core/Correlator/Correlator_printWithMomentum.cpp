#include "Correlator.ih"

namespace Core
{
  template< >
  void Correlator< Dirac::Matrix >::printWithMomentum(std::ostream &out, int const * const momentum) const
  {
    size_t prec_tmp = out.precision();
    // out << "T        real part   imaginary part\n";
    for (size_t t = 0; t < T(); t++)
    {
      out.width(3);
      out << t << "  ";
      out.width(2);
      out << std::showpos << momentum[0] << " ";
      out.width(2);
      out << std::showpos << momentum[1] << " ";
      out.width(2);
      out << std::showpos << momentum[2] << "  ";
      out << std::scientific << std::showpos << std::setprecision(8);
      out << ((*this)[t]).trace().real() << "  " << ((*this)[t]).trace().imag() << std::endl;
    }
    out.precision(prec_tmp);
  }
}