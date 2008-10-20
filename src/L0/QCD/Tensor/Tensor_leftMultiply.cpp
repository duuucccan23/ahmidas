#include "Tensor.ih"

QCD::Tensor &QCD::Tensor::leftMultiply(QCD::Tensor const &other)
{
  std::complex< double > result[12];
  for (size_t colIdx = 0; colIdx < 12; ++colIdx)
  {
    std::fill_n(reinterpret_cast< double* >(result), 24, 0.0);
    for (size_t resIdx = 0; resIdx < 12; ++resIdx)
      for (size_t rowIdx = 0; rowIdx < 12; ++rowIdx)
        result[resIdx] += other.d_data[rowIdx + 12 * resIdx] * d_data[colIdx + 12 * rowIdx];
    for (size_t resIdx = 0; resIdx < 12; ++resIdx)
      d_data[colIdx + 12 * resIdx] = result[resIdx];
  }
}
