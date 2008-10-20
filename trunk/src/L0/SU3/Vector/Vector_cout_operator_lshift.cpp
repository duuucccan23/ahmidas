#include "Vector.ih"

std::ostream &SU3::operator<<(std::ostream &out, SU3::Vector const &vec)
{
  out << std::scientific << std::setprecision(2)
      << "[ " << ((vec.d_data[0].real() > 0) ? "  " : "- ")
              << std::abs(vec.d_data[0].real())
              << ((vec.d_data[0].imag() > 0) ? " + " : " - ") 
              << std::abs(vec.d_data[0].imag()) << " * i    "
              << ((vec.d_data[1].real() > 0) ? "  " : "- ") 
              << std::abs(vec.d_data[1].real())
              << ((vec.d_data[1].imag() > 0) ? " + " : " - ") 
              << std::abs(vec.d_data[1].imag()) << " * i    "
              << ((vec.d_data[2].real() > 0) ? "  " : "- ") 
              << std::abs(vec.d_data[2].real())
              << ((vec.d_data[2].imag() > 0) ? " + " : " - ") 
              << std::abs(vec.d_data[2].imag()) << " * i  ]"
      	    << std::endl;
  return out;
}
