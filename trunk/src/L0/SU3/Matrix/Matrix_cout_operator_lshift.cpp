#include "Matrix.ih"

std::ostream &SU3::operator<<(std::ostream &out, SU3::Matrix const &mat)
{
  out << std::scientific << std::setprecision(3) << std::showpos
      << "[ " << mat.d_data[0].real() <<  "  " << mat.d_data[0].imag() << " * i   "
              << mat.d_data[1].real() <<  "  " << mat.d_data[1].imag() << " * i   "
              << mat.d_data[2].real() <<  "  " << mat.d_data[2].imag() << " * i ]\n"

      << "[ " << mat.d_data[3].real() <<  "  " << mat.d_data[3].imag() << " * i   "
              << mat.d_data[4].real() <<  "  " << mat.d_data[4].imag() << " * i   "
              << mat.d_data[5].real() <<  "  " << mat.d_data[5].imag() << " * i ]\n"

      << "[ " << mat.d_data[6].real() <<  "  " << mat.d_data[6].imag() << " * i   "
              << mat.d_data[7].real() <<  "  " << mat.d_data[7].imag() << " * i   "
              << mat.d_data[8].real() <<  "  " << mat.d_data[8].imag() << " * i ]\n"
 	    << std::endl;
  return out;
}
