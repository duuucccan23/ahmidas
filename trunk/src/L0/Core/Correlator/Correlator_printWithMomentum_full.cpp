#include "Correlator.ih"

namespace Core
{

   // NOTE: this convention seems to differ from an older one by a transpose(), should be fixed shortly
  template< >
  void Correlator< Dirac::Matrix >::printWithMomentum_full(std::ostream &out, int const * const momentum, std::string const& prefix) const
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

      out.precision(8);
      out   << std::scientific << std::showpos
            << ((*this)[t])[ 0].real() << "  " << ((*this)[t])[ 0].imag() << "  "
            << ((*this)[t])[ 1].real() << "  " << ((*this)[t])[ 1].imag() << "  "
            << ((*this)[t])[ 2].real() << "  " << ((*this)[t])[ 2].imag() << "  "
            << ((*this)[t])[ 3].real() << "  " << ((*this)[t])[ 3].imag() << "  "
            << std::endl;

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

      out   << std::scientific << std::showpos
            << ((*this)[t])[ 4].real() << "  " << ((*this)[t])[ 4].imag() << "  "
            << ((*this)[t])[ 5].real() << "  " << ((*this)[t])[ 5].imag() << "  "
            << ((*this)[t])[ 6].real() << "  " << ((*this)[t])[ 6].imag() << "  "
            << ((*this)[t])[ 7].real() << "  " << ((*this)[t])[ 7].imag() << "  "
            << std::endl;

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

      out   << std::scientific << std::showpos
            << ((*this)[t])[ 8].real() << "  " << ((*this)[t])[ 8].imag() << "  "
            << ((*this)[t])[ 9].real() << "  " << ((*this)[t])[ 9].imag() << "  "
            << ((*this)[t])[10].real() << "  " << ((*this)[t])[10].imag() << "  "
            << ((*this)[t])[11].real() << "  " << ((*this)[t])[11].imag() << "  "
            << std::endl;

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

      out   << std::scientific << std::showpos
            << ((*this)[t])[12].real() << "  " << ((*this)[t])[12].imag() << "  "
            << ((*this)[t])[13].real() << "  " << ((*this)[t])[13].imag() << "  "
            << ((*this)[t])[14].real() << "  " << ((*this)[t])[14].imag() << "  "
            << ((*this)[t])[15].real() << "  " << ((*this)[t])[15].imag() << "  "
            << std::endl;
     out.unsetf(std::ios::showpos);
    }
    out.precision(prec_tmp);
  }
}