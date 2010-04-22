#include "Matrix.ih"

namespace Dirac
{
  std::ostream &operator<<(std::ostream &out, Matrix const &mat)
  {
    out << std::scientific << std::setprecision(3) << std::showpos
        << "[ " << mat.d_data[ 0].real() << "  " << mat.d_data[ 0].imag() << " * i   "
                << mat.d_data[ 1].real() << "  " << mat.d_data[ 1].imag() << " * i   "
                << mat.d_data[ 2].real() << "  " << mat.d_data[ 2].imag() << " * i   "
                << mat.d_data[ 3].real() << "  " << mat.d_data[ 3].imag() << " * i  ]"
                << std::endl
        << "[ " << mat.d_data[ 4].real() << "  " << mat.d_data[ 4].imag() << " * i   "
                << mat.d_data[ 5].real() << "  " << mat.d_data[ 5].imag() << " * i   "
                << mat.d_data[ 6].real() << "  " << mat.d_data[ 6].imag() << " * i   "
                << mat.d_data[ 7].real() << "  " << mat.d_data[ 7].imag() << " * i  ]"
                << std::endl
        << "[ " << mat.d_data[ 8].real() << "  " << mat.d_data[ 8].imag() << " * i   "
                << mat.d_data[ 9].real() << "  " << mat.d_data[ 9].imag() << " * i   "
                << mat.d_data[10].real() << "  " << mat.d_data[10].imag() << " * i   "
                << mat.d_data[11].real() << "  " << mat.d_data[11].imag() << " * i  ]"
                << std::endl
        << "[ " << mat.d_data[12].real() << "  " << mat.d_data[12].imag() << " * i   "
                << mat.d_data[13].real() << "  " << mat.d_data[13].imag() << " * i   "
                << mat.d_data[14].real() << "  " << mat.d_data[14].imag() << " * i   "
                << mat.d_data[15].real() << "  " << mat.d_data[15].imag() << " * i  ]"
                << std::endl;
    return out;
  }
}
