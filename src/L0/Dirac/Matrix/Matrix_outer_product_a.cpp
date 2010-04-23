#include "Matrix.ih"

namespace Dirac
{

  void Matrix::outer_product(Matrix const &other, std::complex< double > * const result, OuterProductIndexOrder const idxOrd,
        Base::DiracIndex sourceIndex, Base::DiracIndex sinkIndex) const
  {

  size_t const offset(4*sourceIndex);

    switch (idxOrd)
    {
      case order_OUTER_FIXED:
        // this is the outer product of a row of *this (specified by sinkIndex) and
        // a column of other (specified by sourceIndex)
        std::transform(other.d_data + offset, other.d_data + offset + 4, result,
                      std::bind1st(std::multiplies< std::complex< double > >(), (this->d_data)[sinkIndex     ]));
        std::transform(other.d_data + offset, other.d_data + offset + 4, result +  4,
                      std::bind1st(std::multiplies< std::complex< double > >(), (this->d_data)[sinkIndex +  4]));
        std::transform(other.d_data + offset, other.d_data + offset + 4, result +  8,
                      std::bind1st(std::multiplies< std::complex< double > >(), (this->d_data)[sinkIndex +  8]));
        std::transform(other.d_data + offset, other.d_data + offset + 4, result + 12,
                      std::bind1st(std::multiplies< std::complex< double > >(), (this->d_data)[sinkIndex + 12]));
        break;
      case order_FIRST_FIXED:
        std::transform(other.d_data, other.d_data + 16, result,
                       std::bind1st(std::multiplies< std::complex< double > >(), (this->d_data)[sinkIndex + offset]));
        break;
      case order_FIRST_OUTER_DELTA:
      {
        // this could be done more efficiently, since we only use one column of the matrix product
        Dirac::Matrix tmp1(other);
        tmp1 *= (*this);
        std::fill_n(result, 16, std::complex< double >(0, 0));
        result[sourceIndex     ] = tmp1.d_data[sinkIndex     ];
        result[sourceIndex +  4] = tmp1.d_data[sinkIndex +  4];
        result[sourceIndex +  8] = tmp1.d_data[sinkIndex +  8];
        result[sourceIndex + 12] = tmp1.d_data[sinkIndex + 12];
        break;
      }
      case order_SECOND_OUTER_DELTA:
      {
        Dirac::Matrix tmp2(*this);
        tmp2 *= other;
        std::fill_n(result, 16, std::complex< double >(0, 0));
        std::copy(tmp2.d_data + 4*sinkIndex, tmp2.d_data + 4*sinkIndex + 4, result + offset);
        break;
      }
      case order_BOTH_OUTER_DELTA:
        std::fill_n(result, 16, std::complex< double >(0, 0));
        result[sinkIndex + offset] = ((*this)*other).trace();
        break;
//       default:
//         std::fill_n(result, 16, std::complex< double >(0, 0));
     }
  }
}
