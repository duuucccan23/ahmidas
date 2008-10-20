#include "Tensor.ih"

QCD::Tensor::Tensor(QCD::hcTensor const &other)
{
  for (size_t rowIdx = 0; rowIdx < 12; ++rowIdx)
    for (size_t colIdx = 0; colIdx < 12; ++colIdx)
      d_data[rowIdx + 12 * colIdx] = std::conj(other.d_parent.d_data[colIdx + 12 * rowIdx]);
}
