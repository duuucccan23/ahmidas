#include "Correlator.ih"

namespace Core
{
  template< >
  void Correlator< Dirac::Matrix >::printWithMomentum(std::ostream &out, int const * const momentum, std::string const& prefix) const
  {
    size_t prec_tmp = out.precision();
    // out << "T        real part   imaginary part\n";
    for (size_t t = 0; t < T(); t++)
    {
      out << prefix << " ";
      out << std::setw(3)    << std::noshowpos << t           << "  ";
      out.setf(std::ios::showpos);

      if (momentum[0] == 0)
        out << "+0 ";
      else
        out << std::setw(2)    << std::showpos   << momentum[0] << " ";

      if (momentum[1] == 0)
        out << "+0 ";
      else
        out << std::setw(2)    << std::showpos   << momentum[1] << " ";

      if (momentum[2] == 0)
        out << "+0 ";
      else
        out << std::setw(2)    << std::showpos   << momentum[2] << "  ";

      out << std::scientific << std::showpos   << std::setprecision(8);
      out << ((*this)[t]).trace().real() << "  "
          << ((*this)[t]).trace().imag() << std::endl;
      out.unsetf(std::ios::showpos);
    }
    out.precision(prec_tmp);
  }
}