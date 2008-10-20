#include "Tensor.ih"

QCD::Tensor &QCD::Tensor::rightMultiply(QCD::Tensor const &other)
{
  std::complex< double > result[12];
  for (size_t rowIdx = 0; rowIdx < 12; ++rowIdx)
  {
    std::fill_n(reinterpret_cast< double* >(result), 24, 0.0);
    for (size_t resIdx = 0; resIdx < 12; ++resIdx)
      for (size_t colIdx = 0; colIdx < 12; ++colIdx)
        result[resIdx] += d_data[colIdx + 12 * rowIdx] * other.d_data[resIdx + 12 * colIdx];
    std::copy(reinterpret_cast< double* >(result), reinterpret_cast< double* >(result) + 24,
              reinterpret_cast< double* >(d_data) + 24 * rowIdx);
  }
}

