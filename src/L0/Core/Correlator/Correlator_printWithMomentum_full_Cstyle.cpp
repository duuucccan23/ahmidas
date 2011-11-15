#include "Correlator.ih"
#include <cstdio>
namespace Core
{

  template< >
  void Correlator< Dirac::Matrix >::printWithMomentum_full_Cstyle(FILE * out, int const * const momentum, std::string const& prefix) const
  {
    for (int t = 0; t < int(T()); t++)
    {
      fprintf(out,"%s %+3d %+2d %+2d %+2d  %9.17e %9.17e %9.17e %9.17e %9.17e %9.17e %9.17e %9.17e\n", prefix.c_str(), t,
              momentum[0], momentum[1], momentum[2],
              ((*this)[t])[ 0].real(), ((*this)[t])[ 0].imag(),
              ((*this)[t])[ 1].real(), ((*this)[t])[ 1].imag(),
              ((*this)[t])[ 2].real(), ((*this)[t])[ 2].imag(),
              ((*this)[t])[ 3].real(), ((*this)[t])[ 3].imag()); 
      fprintf(out,"%s %+3d %+2d %+2d %+2d  %9.17e %9.17e %9.17e %9.17e %9.17e %9.17e %9.17e %9.17e\n", prefix.c_str(), t,
              momentum[0], momentum[1], momentum[2],
              ((*this)[t])[ 4].real(), ((*this)[t])[ 4].imag(),
              ((*this)[t])[ 5].real(), ((*this)[t])[ 5].imag(),
              ((*this)[t])[ 6].real(), ((*this)[t])[ 6].imag(),
              ((*this)[t])[ 7].real(), ((*this)[t])[ 7].imag()); 
      fprintf(out,"%s %+3d %+2d %+2d %+2d  %9.17e %9.17e %9.17e %9.17e %9.17e %9.17e %9.17e %9.17e\n", prefix.c_str(), t,
              momentum[0], momentum[1], momentum[2],
             ((*this)[t])[ 8].real(), ((*this)[t])[ 8].imag(),
             ((*this)[t])[ 9].real(), ((*this)[t])[ 9].imag(),
             ((*this)[t])[10].real(), ((*this)[t])[10].imag(),
             ((*this)[t])[11].real(), ((*this)[t])[11].imag()); 
      fprintf(out,"%s %+3d %+2d %+2d %+2d  %9.17e %9.17e %9.17e %9.17e %9.17e %9.17e %9.17e %9.17e\n", prefix.c_str(), t,
              momentum[0], momentum[1], momentum[2],
             ((*this)[t])[12].real(), ((*this)[t])[12].imag(),
             ((*this)[t])[13].real(), ((*this)[t])[13].imag(),
             ((*this)[t])[14].real(), ((*this)[t])[14].imag(),
             ((*this)[t])[15].real(), ((*this)[t])[15].imag()); 

    }
  }
}
