#include "Matrix.ih"

namespace Dirac
{
  void Matrix::operator*=(Matrix const &rhs)
  {
    Matrix tmpR(rhs);
    Matrix tmpL(*this);
    std::complex< double > zero(0, 0);
    //transpose in order to use inner_product of ranges
    std::swap(tmpR.d_data[ 1], tmpR.d_data[ 4]);
    std::swap(tmpR.d_data[ 2], tmpR.d_data[ 8]);
    std::swap(tmpR.d_data[ 3], tmpR.d_data[12]);
    std::swap(tmpR.d_data[ 6], tmpR.d_data[ 9]);
    std::swap(tmpR.d_data[ 7], tmpR.d_data[13]);
    std::swap(tmpR.d_data[11], tmpR.d_data[14]);
    d_data[ 0] = std::inner_product(tmpL.d_data     , tmpL.d_data +  4, tmpR.d_data,      zero);
    d_data[ 1] = std::inner_product(tmpL.d_data     , tmpL.d_data +  4, tmpR.d_data +  4, zero);
    d_data[ 2] = std::inner_product(tmpL.d_data     , tmpL.d_data +  4, tmpR.d_data +  8, zero);
    d_data[ 3] = std::inner_product(tmpL.d_data     , tmpL.d_data +  4, tmpR.d_data + 12, zero);
    d_data[ 4] = std::inner_product(tmpL.d_data +  4, tmpL.d_data +  8, tmpR.d_data,      zero);
    d_data[ 5] = std::inner_product(tmpL.d_data +  4, tmpL.d_data +  8, tmpR.d_data +  4, zero);
    d_data[ 6] = std::inner_product(tmpL.d_data +  4, tmpL.d_data +  8, tmpR.d_data +  8, zero);
    d_data[ 7] = std::inner_product(tmpL.d_data +  4, tmpL.d_data +  8, tmpR.d_data + 12, zero);
    d_data[ 8] = std::inner_product(tmpL.d_data +  8, tmpL.d_data + 12, tmpR.d_data,      zero);
    d_data[ 9] = std::inner_product(tmpL.d_data +  8, tmpL.d_data + 12, tmpR.d_data +  4, zero);
    d_data[10] = std::inner_product(tmpL.d_data +  8, tmpL.d_data + 12, tmpR.d_data +  8, zero);
    d_data[11] = std::inner_product(tmpL.d_data +  8, tmpL.d_data + 12, tmpR.d_data + 12, zero);
    d_data[12] = std::inner_product(tmpL.d_data + 12, tmpL.d_data + 16, tmpR.d_data,      zero);
    d_data[13] = std::inner_product(tmpL.d_data + 12, tmpL.d_data + 16, tmpR.d_data +  4, zero);
    d_data[14] = std::inner_product(tmpL.d_data + 12, tmpL.d_data + 16, tmpR.d_data +  8, zero);
    d_data[15] = std::inner_product(tmpL.d_data + 12, tmpL.d_data + 16, tmpR.d_data + 12, zero);
  }
}