#include "Matrix.ih"

std::ostream &SU3::operator<<(std::ostream &out, SU3::Matrix const &mat)
{
  out << std::scientific << std::setprecision(2)
      << "[ " << ((mat.d_data[0].real() > 0) ? "  " : "- ")
              << std::abs(mat.d_data[0].real())
              << ((mat.d_data[0].imag() > 0) ? " + " : " - ") 
              << std::abs(mat.d_data[0].imag()) << " * i    "
              << ((mat.d_data[1].real() > 0) ? "  " : "- ") 
              << std::abs(mat.d_data[1].real())
              << ((mat.d_data[1].imag() > 0) ? " + " : " - ") 
              << std::abs(mat.d_data[1].imag()) << " * i    "
              << ((mat.d_data[2].real() > 0) ? "  " : "- ") 
              << std::abs(mat.d_data[2].real())
              << ((mat.d_data[2].imag() > 0) ? " + " : " - ") 
              << std::abs(mat.d_data[2].imag()) << " * i  ]\n"
      
      << "[ " << ((mat.d_data[3].real() > 0) ? "  " : "- ")
              << std::abs(mat.d_data[3].real())
              << ((mat.d_data[3].imag() > 0) ? " + " : " - ") 
              << std::abs(mat.d_data[3].imag()) << " * i    "
              << ((mat.d_data[4].real() > 0) ? "  " : "- ") 
              << std::abs(mat.d_data[4].real())
              << ((mat.d_data[4].imag() > 0) ? " + " : " - ") 
              << std::abs(mat.d_data[4].imag()) << " * i    "
              << ((mat.d_data[5].real() > 0) ? "  " : "- ") 
              << std::abs(mat.d_data[5].real())
              << ((mat.d_data[5].imag() > 0) ? " + " : " - ") 
              << std::abs(mat.d_data[5].imag()) << " * i  ]\n"
              
      << "[ " << ((mat.d_data[6].real() > 0) ? "  " : "- ")
              << std::abs(mat.d_data[6].real())
              << ((mat.d_data[6].imag() > 0) ? " + " : " - ") 
              << std::abs(mat.d_data[6].imag()) << " * i    "
              << ((mat.d_data[7].real() > 0) ? "  " : "- ") 
              << std::abs(mat.d_data[7].real())
              << ((mat.d_data[7].imag() > 0) ? " + " : " - ") 
              << std::abs(mat.d_data[7].imag()) << " * i    "
              << ((mat.d_data[8].real() > 0) ? "  " : "- ") 
              << std::abs(mat.d_data[8].real())
              << ((mat.d_data[8].imag() > 0) ? " + " : " - ") 
              << std::abs(mat.d_data[8].imag()) << " * i  ]"
	    << std::endl;
  return out;
}
