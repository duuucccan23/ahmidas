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
      out << std::setw(3)    << std::noshowpos << t           << "  ";
      out.setf(std::ios::showpos);
      out << std::setw(2)    << std::showpos   << momentum[0] << " ";
      out << std::setw(2)    << std::showpos   << momentum[1] << " ";
      out << std::setw(2)    << std::showpos   << momentum[2] << "  ";
      out << std::scientific << std::showpos   << std::setprecision(8);
      out << ((*this)[t]).trace().real() << "  "
          << ((*this)[t]).trace().imag() << std::endl;
      out.unsetf(std::ios::showpos);
    }
    out.precision(prec_tmp);
  }
}