#include "Tensor.ih"

QCD::Tensor &QCD::Tensor::rightMultiply(QCD::hcTensor const &other)
{
  std::complex< double > result[144];
  std::fill_n(result, 144, std::complex< double >(0.0, 0.0)); // Default construction is not sufficient...
  for (size_t rowIdx = 0; rowIdx < 12; ++rowIdx)
    for (size_t colIdx = 0; colIdx < 12; ++colIdx)
      for (size_t resIdx = 0; resIdx < 12; ++resIdx)
        result[rowIdx + colIdx * 12] += d_data[resIdx + 12 * colIdx] * other[rowIdx + 12 * resIdx];
  std::copy(result, result + 144, d_data);
  return *this;
}
