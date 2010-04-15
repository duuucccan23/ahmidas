#include "Matrix.ih"

namespace Dirac
{

//   void Matrix::outer_product(Matrix const &other, Matrix* result) const
//   {
//     std::transform(d_data, d_data + 16, (*(result     )).d_data,
//                    std::bind1st(std::multiplies< std::complex< double > >(), *(other.d_data     )));
//     std::transform(d_data, d_data + 16, (*(result +  1)).d_data,
//                    std::bind1st(std::multiplies< std::complex< double > >(), *(other.d_data +  1)));
//     std::transform(d_data, d_data + 16, (*(result +  2)).d_data,
//                    std::bind1st(std::multiplies< std::complex< double > >(), *(other.d_data +  2)));
//     std::transform(d_data, d_data + 16, (*(result +  3)).d_data,
//                    std::bind1st(std::multiplies< std::complex< double > >(), *(other.d_data +  3)));
//     std::transform(d_data, d_data + 16, (*(result +  4)).d_data,
//                    std::bind1st(std::multiplies< std::complex< double > >(), *(other.d_data +  4)));
//     std::transform(d_data, d_data + 16, (*(result +  5)).d_data,
//                    std::bind1st(std::multiplies< std::complex< double > >(), *(other.d_data +  5)));
//     std::transform(d_data, d_data + 16, (*(result +  6)).d_data,
//                    std::bind1st(std::multiplies< std::complex< double > >(), *(other.d_data +  6)));
//     std::transform(d_data, d_data + 16, (*(result +  7)).d_data,
//                    std::bind1st(std::multiplies< std::complex< double > >(), *(other.d_data +  7)));
//     std::transform(d_data, d_data + 16, (*(result +  8)).d_data,
//                    std::bind1st(std::multiplies< std::complex< double > >(), *(other.d_data +  8)));
//     std::transform(d_data, d_data + 16, (*(result +  9)).d_data,
//                    std::bind1st(std::multiplies< std::complex< double > >(), *(other.d_data +  9)));
//     std::transform(d_data, d_data + 16, (*(result + 10)).d_data,
//                    std::bind1st(std::multiplies< std::complex< double > >(), *(other.d_data + 10)));
//     std::transform(d_data, d_data + 16, (*(result + 11)).d_data,
//                    std::bind1st(std::multiplies< std::complex< double > >(), *(other.d_data + 11)));
//     std::transform(d_data, d_data + 16, (*(result + 12)).d_data,
//                    std::bind1st(std::multiplies< std::complex< double > >(), *(other.d_data + 12)));
//     std::transform(d_data, d_data + 16, (*(result + 13)).d_data,
//                    std::bind1st(std::multiplies< std::complex< double > >(), *(other.d_data + 13)));
//     std::transform(d_data, d_data + 16, (*(result + 14)).d_data,
//                    std::bind1st(std::multiplies< std::complex< double > >(), *(other.d_data + 14)));
//     std::transform(d_data, d_data + 16, (*(result + 15)).d_data,
//                    std::bind1st(std::multiplies< std::complex< double > >(), *(other.d_data + 15)));
//   }


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
        // column 0 of other times  4th element of this->d_data
        std::transform(other.d_data,      other.d_data +  4, (result +   4),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  4)));
        // column 0 of other times  8th element of this->d_data
        std::transform(other.d_data,      other.d_data +  4, (result +   8),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  8)));
        // column 0 of other times 12th element of this->d_data
        std::transform(other.d_data,      other.d_data +  4, (result +  12),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 12)));
        // column 1 of other times  0th element of this->d_data, and so on ...
        std::transform(other.d_data +  4, other.d_data +  8, (result +  16),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data     )));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  20),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  4)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  24),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  8)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  36),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 12)));
        // column 2 of other times  0th element of this->d_data, and so on ...
        std::transform(other.d_data +  8, other.d_data + 12, (result +  40),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data     )));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  44),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  4)));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  48),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  8)));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  52),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 12)));
        // column 3 of other times  0th element of this->d_data, and so on ...
        std::transform(other.d_data + 12, other.d_data + 16, (result +  40),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data     )));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  44),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  4)));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  48),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  8)));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  52),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 12)));
        // column 0 of other times  1st element of this->d_data, and so on ...
        std::transform(other.d_data,      other.d_data +  4, (result      ),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  1)));
        std::transform(other.d_data,      other.d_data +  4, (result +   4),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  5)));
        std::transform(other.d_data,      other.d_data +  4, (result +   8),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  9)));
        std::transform(other.d_data,      other.d_data +  4, (result +  12),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 13)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  16),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  1)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  20),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  5)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  24),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  9)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  36),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 13)));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  40),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  1)));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  44),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  5)));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  48),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  9)));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  52),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 13)));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  40),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  1)));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  44),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  5)));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  48),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  9)));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  52),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 13)));
        // column 0 of other times 2nd element of this->d_data, and so on ...
        std::transform(other.d_data,      other.d_data +  4, (result      ),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  2)));
        std::transform(other.d_data,      other.d_data +  4, (result +   4),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  6)));
        std::transform(other.d_data,      other.d_data +  4, (result +   8),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 10)));
        std::transform(other.d_data,      other.d_data +  4, (result +  12),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 14)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  16),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  2)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  20),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  6)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  24),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 10)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  36),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 14)));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  40),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  2)));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  44),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  6)));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  48),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 10)));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  52),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 14)));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  40),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  2)));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  44),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  6)));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  48),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 10)));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  52),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 14)));
        // column 0 of other times 3rd element of this->d_data, and so on ...
        std::transform(other.d_data,      other.d_data +  4, (result      ),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  3)));
        std::transform(other.d_data,      other.d_data +  4, (result +   4),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  7)));
        std::transform(other.d_data,      other.d_data +  4, (result +   8),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 11)));
        std::transform(other.d_data,      other.d_data +  4, (result +  12),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 15)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  16),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  3)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  20),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  7)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  24),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 11)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  36),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 15)));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  40),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  3)));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  44),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  7)));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  48),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 11)));
        std::transform(other.d_data +  8, other.d_data + 12, (result +  52),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 15)));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  40),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  3)));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  44),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  7)));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  48),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 11)));
        std::transform(other.d_data + 12, other.d_data + 16, (result +  52),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 15)));
        (const_cast< Dirac::Matrix * >(&other))->transpose();
        break;
      }
  }

}