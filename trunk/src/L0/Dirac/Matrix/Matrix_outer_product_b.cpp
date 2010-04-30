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
        // have to multiply row of other by column of this
        // row no. 0 times column no. 0:
        std::transform(other.d_data,      other.d_data +  4, (result      ),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data     )));
        std::transform(other.d_data,      other.d_data +  4, (result +   4),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  4)));
        std::transform(other.d_data,      other.d_data +  4, (result +   8),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  8)));
        std::transform(other.d_data,      other.d_data +  4, (result +  12),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 12)));
        // row no. 0 times column no. 1:
        std::transform(other.d_data,      other.d_data +  4, (result +  16),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  1)));
        std::transform(other.d_data,      other.d_data +  4, (result +  20),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  5)));
        std::transform(other.d_data,      other.d_data +  4, (result +  24),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  9)));
        std::transform(other.d_data,      other.d_data +  4, (result +  28),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 13)));
        // row no. 0 times column no. 2:
        std::transform(other.d_data,      other.d_data +  4, (result +  32),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  2)));
        std::transform(other.d_data,      other.d_data +  4, (result +  36),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  6)));
        std::transform(other.d_data,      other.d_data +  4, (result +  40),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 10)));
        std::transform(other.d_data,      other.d_data +  4, (result +  44),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 14)));
        // row no. 0 times column no. 3:
        std::transform(other.d_data,      other.d_data +  4, (result +  48),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  3)));
        std::transform(other.d_data,      other.d_data +  4, (result +  52),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  7)));
        std::transform(other.d_data,      other.d_data +  4, (result +  56),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 11)));
        std::transform(other.d_data,      other.d_data +  4, (result +  60),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 15)));
        // row no. 1 times column no. 0:
        std::transform(other.d_data +  4, other.d_data +  8, (result +  64),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data     )));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  68),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  4)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  72),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  8)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  76),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 12)));
        // row no. 1 times column no. 1:
        std::transform(other.d_data +  4, other.d_data +  8, (result +  80),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  1)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  84),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  5)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  88),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  9)));
        std::transform(other.d_data +  4, other.d_data +  8, (result +  92),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 13)));
        // row no. 1 times column no. 2:
        std::transform(other.d_data +  4, other.d_data +  8, (result +  96),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  2)));
        std::transform(other.d_data +  4, other.d_data +  8, (result + 100),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  6)));
        std::transform(other.d_data +  4, other.d_data +  8, (result + 104),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 10)));
        std::transform(other.d_data +  4, other.d_data +  8, (result + 108),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 14)));
        // row no. 1 times column no. 3:
        std::transform(other.d_data +  4, other.d_data +  8, (result + 112),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  3)));
        std::transform(other.d_data +  4, other.d_data +  8, (result + 116),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  7)));
        std::transform(other.d_data +  4, other.d_data +  8, (result + 120),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 11)));
        std::transform(other.d_data +  4, other.d_data +  8, (result + 124),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 15)));
        // row no. 2 times column no. 0:
        std::transform(other.d_data +  8, other.d_data + 12, (result + 128),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data     )));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 132),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  4)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 136),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  8)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 140),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 12)));
        // row no. 2 times column no. 1:
        std::transform(other.d_data +  8, other.d_data + 12, (result + 144),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  1)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 148),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  5)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 152),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  9)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 156),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 13)));
        // row no. 2 times column no. 2:
        std::transform(other.d_data +  8, other.d_data + 12, (result + 160),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  2)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 164),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  6)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 168),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 10)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 172),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 14)));
        // row no. 2 times column no. 3
        std::transform(other.d_data +  8, other.d_data + 12, (result + 176),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  3)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 180),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  7)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 184),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 11)));
        std::transform(other.d_data +  8, other.d_data + 12, (result + 188),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 15)));
        // row no. 3 times column no. 0
        std::transform(other.d_data + 12, other.d_data + 16, (result + 192),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data     )));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 196),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  4)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 200),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  8)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 204),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 12)));
        // row no. 3 times column no. 1
        std::transform(other.d_data + 12, other.d_data + 16, (result + 208),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  1)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 212),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  5)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 216),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  9)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 220),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 13)));
        // row no. 3 times column no. 2
        std::transform(other.d_data + 12, other.d_data + 16, (result + 224),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  2)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 228),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  6)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 232),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 10)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 236),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 14)));
        // row no. 3 times column no. 3
        std::transform(other.d_data + 12, other.d_data + 16, (result + 240),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  3)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 244),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data +  7)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 248),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 11)));
        std::transform(other.d_data + 12, other.d_data + 16, (result + 252),
                      std::bind1st(std::multiplies< std::complex< double > >(), *(d_data + 15)));
        break;
        case order_FIRST_OUTER_DELTA:
        {
          // not tested yet...
          assert(false);
          std::fill_n(result, 256, std::complex< double >(0, 0));
          Dirac::Matrix const tmp((*this)*other);
          std::copy(tmp.d_data,      tmp.d_data +  4, result);
          std::copy(tmp.d_data,      tmp.d_data +  4, result +  20);
          std::copy(tmp.d_data,      tmp.d_data +  4, result +  40);
          std::copy(tmp.d_data,      tmp.d_data +  4, result +  60);
          std::copy(tmp.d_data +  4, tmp.d_data +  4, result +  65);
          std::copy(tmp.d_data +  4, tmp.d_data +  8, result +  85);
          std::copy(tmp.d_data +  4, tmp.d_data +  8, result + 105);
          std::copy(tmp.d_data +  4, tmp.d_data +  8, result + 125);
          std::copy(tmp.d_data +  8, tmp.d_data +  8, result + 130);
          std::copy(tmp.d_data +  8, tmp.d_data + 12, result + 150);
          std::copy(tmp.d_data +  8, tmp.d_data + 12, result + 170);
          std::copy(tmp.d_data +  8, tmp.d_data + 12, result + 190);
          std::copy(tmp.d_data + 12, tmp.d_data + 16, result + 195);
          std::copy(tmp.d_data + 12, tmp.d_data + 16, result + 215);
          std::copy(tmp.d_data + 12, tmp.d_data + 16, result + 235);
          std::copy(tmp.d_data + 12, tmp.d_data + 16, result + 255);
          break;
        }
        case order_SECOND_OUTER_DELTA:
        {
          // not completed yet...
          assert(false);
          std::fill_n(result, 256, std::complex< double >(0, 0));
          Dirac::Matrix const tmp((*this)*other);
          std::fill_n(result, 16, std::complex< double >(0, 0));
          result[  0] = tmp.d_data[  0];
          result[  4] = tmp.d_data[  4];
          result[  8] = tmp.d_data[  8];
          result[ 12] = tmp.d_data[ 12];
          result[ 16] = tmp.d_data[  1];
          result[ 20] = tmp.d_data[  5];
          result[ 24] = tmp.d_data[  9];
          result[ 28] = tmp.d_data[ 13];
          result[ 32] = tmp.d_data[  2];
          result[ 36] = tmp.d_data[  6];
          result[ 40] = tmp.d_data[ 10];
          result[ 44] = tmp.d_data[ 14];
          result[ 48] = tmp.d_data[  3];
          result[ 52] = tmp.d_data[  7];
          result[ 56] = tmp.d_data[ 11];
          result[ 60] = tmp.d_data[ 15];
          result[ 65] = tmp.d_data[  0];
          result[ 69] = tmp.d_data[  4];
          result[ 73] = tmp.d_data[  8];
          result[ 77] = tmp.d_data[ 12];

          // and so on ...

          break;
        }
        case order_BOTH_OUTER_DELTA:
        {
          // not tested yet...
          assert(false);
          std::fill_n(result, 256, std::complex< double >(0, 0));
          Dirac::Matrix const tmp((*this)*other);
          result[  0] = tmp[ 0];
          result[ 17] = tmp[ 4];
          result[ 34] = tmp[ 8];
          result[ 51] = tmp[12];
          result[ 68] = tmp[ 1];
          result[ 85] = tmp[ 5];
          result[102] = tmp[ 9];
          result[119] = tmp[13];
          result[136] = tmp[ 2];
          result[153] = tmp[ 6];
          result[170] = tmp[10];
          result[187] = tmp[14];
          result[204] = tmp[ 3];
          result[221] = tmp[ 7];
          result[238] = tmp[11];
          result[255] = tmp[15];
          break;
        }

      }
  }

}
