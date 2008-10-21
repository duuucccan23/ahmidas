#include "Tensor.ih"

QCD::Tensor &QCD::Tensor::operator=(QCD::Tensor const &other)
{
  if (this != &other)
    std::copy(other.d_data, other.d_data + 144, d_data);
  return *this;
}
