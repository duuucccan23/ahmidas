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
          std::fill_n(result, 256, std::complex< double >(0, 0));
          Dirac::Matrix const tmp(other*(*this));
          std::copy(tmp.d_data,      tmp.d_data +  4, result);
          std::copy(tmp.d_data,      tmp.d_data +  4, result +  20);
          std::copy(tmp.d_data,      tmp.d_data +  4, result +  40);
          std::copy(tmp.d_data,      tmp.d_data +  4, result +  60);
          std::copy(tmp.d_data +  4, tmp.d_data +  8, result +  64);
          std::copy(tmp.d_data +  4, tmp.d_data +  8, result +  84);
          std::copy(tmp.d_data +  4, tmp.d_data +  8, result + 104);
          std::copy(tmp.d_data +  4, tmp.d_data +  8, result + 124);
          std::copy(tmp.d_data +  8, tmp.d_data + 12, result + 128);
          std::copy(tmp.d_data +  8, tmp.d_data + 12, result + 148);
          std::copy(tmp.d_data +  8, tmp.d_data + 12, result + 168);
          std::copy(tmp.d_data +  8, tmp.d_data + 12, result + 188);
          std::copy(tmp.d_data + 12, tmp.d_data + 16, result + 192);
          std::copy(tmp.d_data + 12, tmp.d_data + 16, result + 212);
          std::copy(tmp.d_data + 12, tmp.d_data + 16, result + 232);
          std::copy(tmp.d_data + 12, tmp.d_data + 16, result + 252);
          break;
        }
        case order_SECOND_OUTER_DELTA:
        {
          std::fill_n(result, 256, std::complex< double >(0, 0));
          Dirac::Matrix const tmp(other*(*this));
          // column no. 0 (copied to column no. 0 of result)
          result[  0] = tmp.d_data[  0];
          result[  4] = tmp.d_data[  4];
          result[  8] = tmp.d_data[  8];
          result[ 12] = tmp.d_data[ 12];
          // column no. 1 (copied to column no. 0 of result)
          result[ 16] = tmp.d_data[  1];
          result[ 20] = tmp.d_data[  5];
          result[ 24] = tmp.d_data[  9];
          result[ 28] = tmp.d_data[ 13];
          // column no. 2 (copied to column no. 0 of result)
          result[ 32] = tmp.d_data[  2];
          result[ 36] = tmp.d_data[  6];
          result[ 40] = tmp.d_data[ 10];
          result[ 44] = tmp.d_data[ 14];
          // column no. 3 (copied to column no. 0 of result)
          result[ 48] = tmp.d_data[  3];
          result[ 52] = tmp.d_data[  7];
          result[ 56] = tmp.d_data[ 11];
          result[ 60] = tmp.d_data[ 15];
          // column no. 0 (copied to column no. 1 of result)
          result[ 65] = tmp.d_data[  0];
          result[ 69] = tmp.d_data[  4];
          result[ 73] = tmp.d_data[  8];
          result[ 77] = tmp.d_data[ 12];
          // column no. 1 (copied to column no. 1 of result)
          result[ 81] = tmp.d_data[  1];
          result[ 85] = tmp.d_data[  5];
          result[ 89] = tmp.d_data[  9];
          result[ 93] = tmp.d_data[ 13];
          // column no. 2 (copied to column no. 1 of result)
          result[ 97] = tmp.d_data[  2];
          result[101] = tmp.d_data[  6];
          result[105] = tmp.d_data[ 10];
          result[109] = tmp.d_data[ 14];
          // column no. 3 (copied to column no. 1 of result)
          result[113] = tmp.d_data[  3];
          result[117] = tmp.d_data[  7];
          result[121] = tmp.d_data[ 11];
          result[125] = tmp.d_data[ 15];
          // column no. 0 (copied to column no. 2 of result)
          result[130] = tmp.d_data[  0];
          result[134] = tmp.d_data[  4];
          result[138] = tmp.d_data[  8];
          result[142] = tmp.d_data[ 12];
          // column no. 1 (copied to column no. 2 of result)
          result[146] = tmp.d_data[  1];
          result[150] = tmp.d_data[  5];
          result[154] = tmp.d_data[  9];
          result[158] = tmp.d_data[ 13];
          // column no. 2 (copied to column no. 2 of result)
          result[162] = tmp.d_data[  2];
          result[166] = tmp.d_data[  6];
          result[170] = tmp.d_data[ 10];
          result[174] = tmp.d_data[ 14];
          // column no. 3 (copied to column no. 2 of result)
          result[178] = tmp.d_data[  3];
          result[182] = tmp.d_data[  7];
          result[186] = tmp.d_data[ 11];
          result[190] = tmp.d_data[ 15];
          // column no. 0 (copied to column no. 3 of result)
          result[195] = tmp.d_data[  0];
          result[199] = tmp.d_data[  4];
          result[203] = tmp.d_data[  8];
          result[207] = tmp.d_data[ 12];
          // column no. 1 (copied to column no. 3 of result)
          result[211] = tmp.d_data[  1];
          result[215] = tmp.d_data[  5];
          result[219] = tmp.d_data[  9];
          result[223] = tmp.d_data[ 13];
          // column no. 2 (copied to column no. 3 of result)
          result[227] = tmp.d_data[  2];
          result[231] = tmp.d_data[  6];
          result[235] = tmp.d_data[ 10];
          result[239] = tmp.d_data[ 14];
          // column no. 3 (copied to column no. 3 of result)
          result[243] = tmp.d_data[  3];
          result[247] = tmp.d_data[  7];
          result[251] = tmp.d_data[ 11];
          result[255] = tmp.d_data[ 15];
          break;
        }
        case order_BOTH_OUTER_DELTA:
        {
          std::fill_n(result, 256, std::complex< double >(0, 0));
          std::complex< double > const tmp(((*this)*other).trace());
          result[  0] = tmp;
          result[ 20] = tmp;
          result[ 40] = tmp;
          result[ 60] = tmp;
          result[ 65] = tmp;
          result[ 85] = tmp;
          result[105] = tmp;
          result[125] = tmp;
          result[130] = tmp;
          result[150] = tmp;
          result[170] = tmp;
          result[190] = tmp;
          result[195] = tmp;
          result[215] = tmp;
          result[235] = tmp;
          result[255] = tmp;
          break;
        }

      }
  }

}
