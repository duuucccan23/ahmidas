#include "Tensor.ih"

namespace QCD
{

  void Tensor::getDiracMatrix(Dirac::Matrix &dMatrix, Base::ColourIndex const colour_src, Base::ColourIndex const colour_snk) const
  {
    size_t index = 12*colour_src + colour_snk;
    dMatrix[ 0] = d_data[index      ];
    dMatrix[ 1] = d_data[index +   3];
    dMatrix[ 2] = d_data[index +   6];
    dMatrix[ 3] = d_data[index +   9];
    dMatrix[ 4] = d_data[index +  36];
    dMatrix[ 5] = d_data[index +  39];
    dMatrix[ 6] = d_data[index +  42];
    dMatrix[ 7] = d_data[index +  45];
    dMatrix[ 8] = d_data[index +  72];
    dMatrix[ 9] = d_data[index +  75];
    dMatrix[10] = d_data[index +  78];
    dMatrix[11] = d_data[index +  81];
    dMatrix[12] = d_data[index + 108];
    dMatrix[13] = d_data[index + 111];
    dMatrix[14] = d_data[index + 114];
    dMatrix[15] = d_data[index + 117];
  }
}
