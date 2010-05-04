#include "Matrix.ih"

namespace Dirac
{

  void Matrix::outer_product(Matrix const &other, std::complex< double > * const result, OuterProductIndexOrder const idxOrd,
        Base::DiracIndex sourceIndex, Base::DiracIndex sinkIndex) const
  {

    size_t const offsetSrc(4*sourceIndex);

    switch (idxOrd)
    {
      case order_OUTER_FIXED:
      {
        // this is the outer product of a row of *this (specified by sourceIndex) and
        // a column of other (specified by sinkIndex)
        std::transform(other.d_data + offsetSrc, other.d_data + offsetSrc + 4, result,
                      std::bind1st(std::multiplies< std::complex< double > >(), (this->d_data)[sinkIndex     ]));
        std::transform(other.d_data + offsetSrc, other.d_data + offsetSrc + 4, result +  4,
                      std::bind1st(std::multiplies< std::complex< double > >(), (this->d_data)[sinkIndex +  4]));
        std::transform(other.d_data + offsetSrc, other.d_data + offsetSrc + 4, result +  8,
                      std::bind1st(std::multiplies< std::complex< double > >(), (this->d_data)[sinkIndex +  8]));
        std::transform(other.d_data + offsetSrc, other.d_data + offsetSrc + 4, result + 12,
                      std::bind1st(std::multiplies< std::complex< double > >(), (this->d_data)[sinkIndex + 12]));
        break;
      }
      case order_FIRST_FIXED:
        std::transform(other.d_data, other.d_data + 16, result,
                       std::bind1st(std::multiplies< std::complex< double > >(), (this->d_data)[sinkIndex + offsetSrc]));
        break;
      case order_FIRST_OUTER_DELTA:
      {
        // this could be done more efficiently, since we only use one row of the matrix product
        Dirac::Matrix tmp(other);
        tmp *= (*this);
        std::fill_n(result, 16, std::complex< double >(0, 0));
        std::copy(tmp.d_data + offsetSrc, tmp.d_data + offsetSrc + 4, result + 4*sinkIndex);
        break;
      }
      case order_SECOND_OUTER_DELTA:
      {
        // this could be done more efficiently, since we only use one column of the matrix product
        Dirac::Matrix tmp(other);
        tmp *= (*this);
        std::fill_n(result, 16, std::complex< double >(0, 0));
        result[sourceIndex     ] = tmp.d_data[sinkIndex     ];
        result[sourceIndex +  4] = tmp.d_data[sinkIndex +  4];
        result[sourceIndex +  8] = tmp.d_data[sinkIndex +  8];
        result[sourceIndex + 12] = tmp.d_data[sinkIndex + 12];
        break;
      }
      case order_BOTH_OUTER_DELTA:
        std::fill_n(result, 16, std::complex< double >(0, 0));
        result[4*sinkIndex + sourceIndex] = ((*this)*other).trace();
        break;
     }
  }
}
