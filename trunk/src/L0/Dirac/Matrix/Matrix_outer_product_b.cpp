#include "Matrix.ih"

namespace Dirac
{

  void Matrix::outer_product(Matrix const &other, std::complex< double > * const result, OuterProductIndexOrder const idxOrd) const
  {
    switch (idxOrd)
    {
      case order_FIRST_FIXED:
        std::transform(other.d_data, other.d_data + 16, (result     ),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data     )));
        std::transform(other.d_data, other.d_data + 16, (result +  16),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  1)));
        std::transform(other.d_data, other.d_data + 16, (result +  32),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  2)));
        std::transform(other.d_data, other.d_data + 16, (result +  48),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  3)));
        std::transform(other.d_data, other.d_data + 16, (result +  64),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  4)));
        std::transform(other.d_data, other.d_data + 16, (result +  80),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  5)));
        std::transform(other.d_data, other.d_data + 16, (result +  96),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  6)));
        std::transform(other.d_data, other.d_data + 16, (result + 112),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  7)));
        std::transform(other.d_data, other.d_data + 16, (result + 128),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  8)));
        std::transform(other.d_data, other.d_data + 16, (result + 144),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  9)));
        std::transform(other.d_data, other.d_data + 16, (result + 160),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 10)));
        std::transform(other.d_data, other.d_data + 16, (result + 176),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 11)));
        std::transform(other.d_data, other.d_data + 16, (result + 192),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 12)));
        std::transform(other.d_data, other.d_data + 16, (result + 208),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 13)));
        std::transform(other.d_data, other.d_data + 16, (result + 224),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 14)));
        std::transform(other.d_data, other.d_data + 16, (result + 240),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 15)));
        break;

      case order_OUTER_FIXED:
        // change colums <--> rows for better performance
        (const_cast< Dirac::Matrix * >(&other))->transpose();
        // column 0 of other times  0th element of this->d_data
        std::transform(other.d_data,      other.d_data +  4, (result      ),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data     )));
        // column 0 of other times  1st element of this->d_data
        std::transform(other.d_data,      other.d_data +  4, (result +   4),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  1)));
        // column 0 of other times  2nd element of this->d_data
        std::transform(other.d_data,      other.d_data +  4, (result +   8),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  2)));
        // column 0 of other times 3rd element of this->d_data
        std::transform(other.d_data,      other.d_data +  4, (result +  12),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  3)));
        // column 1 of other times  0th element of this->d_data, and so on ...
        std::transform(other.d_data +  4, other.d_data +  8, (result +  16),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data     )));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  20),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  1)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  24),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  2)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  28),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  3)));
        // column 2 of other times  0th element of this->d_data, and so on ...
        std::transform(other.d_data +  8, other.d_data + 12, (result +  32),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data     )));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  36),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  1)));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  40),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  2)));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  44),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  3)));
        // column 3 of other times  0th element of this->d_data, and so on ...
        std::transform(other.d_data + 12, other.d_data + 16, (result +  48),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data     )));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  52),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  1)));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  56),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  2)));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  60),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  3)));
        // column 0 of other times  4th element of this->d_data, and so on ...
        std::transform(other.d_data,      other.d_data +  4, (result +  64),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  4)));
        std::transform(other.d_data,      other.d_data +  4, (result +  68),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  5)));
        std::transform(other.d_data,      other.d_data +  4, (result +  72),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  6)));
        std::transform(other.d_data,      other.d_data +  4, (result +  76),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  7)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  80),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  4)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  84),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  5)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  88),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  6)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  92),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  7)));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  96),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  4)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 100),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  5)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 104),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  6)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 108),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  7)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 112),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  4)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 116),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  5)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 120),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  6)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 124),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  7)));
        // column 0 of other times 8th element of this->d_data, and so on ...
        std::transform(other.d_data,      other.d_data +  4, (result + 128),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  8)));
        std::transform(other.d_data,      other.d_data +  4, (result + 132),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  9)));
        std::transform(other.d_data,      other.d_data +  4, (result + 136),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 10)));
        std::transform(other.d_data,      other.d_data +  4, (result + 140),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 11)));
        std::transform(other.d_data +  4, other.d_data +  8, (result + 144),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  8)));
        std::transform(other.d_data +  4, other.d_data +  8, (result + 148),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  9)));
        std::transform(other.d_data +  4, other.d_data +  8, (result + 152),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 10)));
        std::transform(other.d_data +  4, other.d_data +  8, (result + 156),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 11)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 160),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  8)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 164),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  9)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 168),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 10)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 172),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 11)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 176),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  8)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 180),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  9)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 184),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 10)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 188),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 11)));
        // column 0 of other times 12th element of this->d_data, and so on ...
        std::transform(other.d_data,      other.d_data +  4, (result + 192),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 12)));
        std::transform(other.d_data,      other.d_data +  4, (result + 196),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 13)));
        std::transform(other.d_data,      other.d_data +  4, (result + 200),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 14)));
        std::transform(other.d_data,      other.d_data +  4, (result + 204),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 15)));
        std::transform(other.d_data +  4, other.d_data +  8, (result + 208),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 12)));
        std::transform(other.d_data +  4, other.d_data +  8, (result + 212),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 13)));
        std::transform(other.d_data +  4, other.d_data +  8, (result + 216),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 14)));
        std::transform(other.d_data +  4, other.d_data +  8, (result + 220),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 15)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 224),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 12)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 228),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 13)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 232),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 14)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 236),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 15)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 240),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 12)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 244),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 13)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 248),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 14)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 252),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 15)));
        (const_cast< Dirac::Matrix * >(&other))->transpose();
        break;
        case order_FIRST_OUTER_DELTA:
        case order_SECOND_OUTER_DELTA:
        case order_BOTH_OUTER_DELTA:
          std::cerr << "This should not happen in Matrix::outer_product(Matrix &, std::complex< double > *, OuterProductIndexOrder)";
          std::cerr << std::endl;
          exit(1);
      }
  }

}
