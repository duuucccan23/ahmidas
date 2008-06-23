#include "Matrix.ih"

std::ostream &SU3::operator<<(std::ostream &out, SU3::Matrix const &mat)
{
  out << mat.d_data[0] << "  " << mat.d_data [1] << "  " << mat.d_data[2] << "\n";
  out << mat.d_data[3] << "  " << mat.d_data [4] << "  " << mat.d_data[5] << "\n";
  out << mat.d_data[6] << "  " << mat.d_data [7] << "  " << mat.d_data[8] << std::endl;
  return out;
}
