#include "Tensor.ih"

namespace QCD
{

  // this function is an experimental one...
  // right is the pointer to a Tensor array of 16 elements
  void multiplyOuterDiracIndices(Tensor const &left, Tensor const * const right, Tensor * const result)
  {
    (const_cast < QCD::Tensor & > (left)).transposeFull();
    std::complex< double > const ZERO(0, 0);
    for (size_t iDirac=0; iDirac<16; iDirac++)
      std::fill_n(result[iDirac].d_data, 144, ZERO);

    //std::complex< double > const ONE(1, 0);
    //std::fill_n((const_cast < QCD::Tensor & > (left)).d_data, 144, ONE);

    std::complex< double > tmp[12];
    for (size_t iD1=0; iD1<16; iD1+=4) // outer source Dirac index of right
    {
      for (size_t iCl=0; iCl<3; iCl++) // sink colour index of left
      {
        for (size_t iDl=0; iDl<4; iDl++) // sink Dirac index of left
        {
          for (size_t iCr=0; iCr<3; iCr++) // source colour index of right
          {
            for (size_t iDr=0; iDr<4; iDr++) // inner source Dirac index of right
            {
              size_t const resIndex((iDr*3+iCr)*12+iCl);
              assert(resIndex<144);
              size_t const rightIndex((iDr*3+iCr)*12);
              assert(rightIndex<144);
              size_t const leftIndex((iDl*3+iCl)*12);
              assert(leftIndex<144);
              // inner sink Dirac index of right = 0
              std::copy(right[iD1    ].d_data+rightIndex,   right[iD1    ].d_data+rightIndex+3,  tmp);
              std::copy(right[iD1 + 1].d_data+rightIndex,   right[iD1 + 1].d_data+rightIndex+3,  tmp);
              std::copy(right[iD1 + 2].d_data+rightIndex,   right[iD1 + 2].d_data+rightIndex+3,  tmp);
              std::copy(right[iD1 + 3].d_data+rightIndex,   right[iD1 + 3].d_data+rightIndex+3,  tmp);
              (result[iDl+iD1]).d_data[resIndex  ] = std::inner_product(left.d_data + leftIndex, left.d_data + leftIndex + 12, tmp, ZERO);
              // inner sink Dirac index of right = 1
              std::copy(right[iD1    ].d_data+rightIndex+3, right[iD1    ].d_data+rightIndex+6,  tmp);
              std::copy(right[iD1 + 1].d_data+rightIndex+3, right[iD1 + 1].d_data+rightIndex+6,  tmp);
              std::copy(right[iD1 + 2].d_data+rightIndex+3, right[iD1 + 2].d_data+rightIndex+6,  tmp);
              std::copy(right[iD1 + 3].d_data+rightIndex+3, right[iD1 + 3].d_data+rightIndex+6,  tmp);
              (result[iDl+iD1]).d_data[resIndex+3] = std::inner_product(left.d_data + leftIndex, left.d_data + leftIndex + 12, tmp, ZERO);
              // inner sink Dirac index of right = 2
              std::copy(right[iD1    ].d_data+rightIndex+6, right[iD1    ].d_data+rightIndex+9,  tmp);
              std::copy(right[iD1 + 1].d_data+rightIndex+6, right[iD1 + 1].d_data+rightIndex+9,  tmp);
              std::copy(right[iD1 + 2].d_data+rightIndex+6, right[iD1 + 2].d_data+rightIndex+9,  tmp);
              std::copy(right[iD1 + 3].d_data+rightIndex+6, right[iD1 + 3].d_data+rightIndex+9,  tmp);
              (result[iDl+iD1]).d_data[resIndex+6] = std::inner_product(left.d_data + leftIndex, left.d_data + leftIndex + 12, tmp, ZERO);
              // inner sink Dirac index of right = 3
              std::copy(right[iD1    ].d_data+rightIndex+9, right[iD1    ].d_data+rightIndex+12, tmp);
              std::copy(right[iD1 + 1].d_data+rightIndex+9, right[iD1 + 1].d_data+rightIndex+12, tmp);
              std::copy(right[iD1 + 2].d_data+rightIndex+9, right[iD1 + 2].d_data+rightIndex+12, tmp);
              std::copy(right[iD1 + 3].d_data+rightIndex+9, right[iD1 + 3].d_data+rightIndex+12, tmp);
              (result[iDl+iD1]).d_data[resIndex+9] = std::inner_product(left.d_data + leftIndex, left.d_data + leftIndex + 12, tmp, ZERO);
            }
          }
        }
      }
    }
    (const_cast < QCD::Tensor & > (left)).transposeFull();

//     for (size_t iDirac=0; iDirac<16; iDirac++)
//       std::cout << iDirac << ":\n" << result[iDirac] << std::endl;
//     exit(2);
  }
}
