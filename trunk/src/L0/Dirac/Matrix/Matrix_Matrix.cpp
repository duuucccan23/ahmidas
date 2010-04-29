#include "Matrix.ih"

namespace Dirac
{
  Matrix::Matrix(std::complex< double > const * const data)
  {
    std::copy(data, data+16, d_data);
  }
}
