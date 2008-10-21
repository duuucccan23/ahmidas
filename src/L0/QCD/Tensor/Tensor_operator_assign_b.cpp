#include "Tensor.ih"

QCD::Tensor &QCD::Tensor::operator=(QCD::hcTensor const &other)
{
  if (this != &other.d_parent)
  {
    for (size_t idx = 0; idx < 144; ++idx)
      d_data[idx] = other[idx];
    return *this;
  }

  std::complex< double > data[144];
  for (size_t idx = 0; idx < 144; ++idx)
    data[idx] = other[idx];
  std::copy(data, data + 144, d_data);
  return *this;
}
