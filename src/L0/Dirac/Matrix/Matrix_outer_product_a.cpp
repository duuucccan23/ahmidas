#include "Matrix.ih"

namespace Dirac
{

  void Matrix::outer_product(Matrix const &other, std::complex< double > * const result, OuterProductIndexOrder const idxOrd,
        Base::DiracIndex sourceIndex, Base::DiracIndex sinkIndex) const
  {

    size_t const offsetSrc(4*sourceIndex);
    size_t const offsetSnk(4*sinkIndex);

    switch (idxOrd)
    {
      case order_OUTER_FIXED:
      (const_cast< Dirac::Matrix * >(&other))->transpose();
        // this is the outer product of a row of *this (specified by sourceIndex) and
        // a column of other (specified by sinkIndex)
        std::transform(other.d_data + offsetSnk, other.d_data + offsetSnk + 4, result,
                      std::bind1st(std::multiplies< std::complex< double > >(), (this->d_data)[offsetSrc    ]));
        std::transform(other.d_data + offsetSnk, other.d_data + offsetSnk + 4, result +  4,
                      std::bind1st(std::multiplies< std::complex< double > >(), (this->d_data)[offsetSrc + 1]));
        std::transform(other.d_data + offsetSnk, other.d_data + offsetSnk + 4, result +  8,
                      std::bind1st(std::multiplies< std::complex< double > >(), (this->d_data)[offsetSrc + 2]));
        std::transform(other.d_data + offsetSnk, other.d_data + offsetSnk + 4, result + 12,
                      std::bind1st(std::multiplies< std::complex< double > >(), (this->d_data)[offsetSrc + 3]));
        (const_cast< Dirac::Matrix * >(&other))->transpose();
        break;
      case order_FIRST_FIXED:
        std::transform(other.d_data, other.d_data + 16, result,
                       std::bind1st(std::multiplies< std::complex< double > >(), (this->d_data)[sinkIndex + offsetSrc]));
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
        std::copy(tmp2.d_data + 4*sinkIndex, tmp2.d_data + 4*sinkIndex + 4, result + offsetSrc);
        break;
      }
      case order_BOTH_OUTER_DELTA:
        std::fill_n(result, 16, std::complex< double >(0, 0));
        result[sinkIndex + offsetSrc] = ((*this)*other).trace();
        break;
//       default:
//         std::fill_n(result, 16, std::complex< double >(0, 0));
     }
  }
}
