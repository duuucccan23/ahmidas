#include "Tensor.ih"

QCD::Tensor &QCD::Tensor::leftMultiply(SU3::Matrix const &mat)
{
  std::complex< double > result0[144];
  std::complex< double > result1[144];
  std::complex< double > result2[144];

  for (size_t idx=0; idx<144; idx+=36)
  {
    std::transform(d_data + idx,      d_data + idx + 12, result0 + idx,
                   std::bind1st(std::multiplies< std::complex< double > >(), mat.d_data[0]));
    std::transform(d_data + idx + 12, d_data + idx + 24, result1 + idx,
                   std::bind1st(std::multiplies< std::complex< double > >(), mat.d_data[1]));
    std::transform(d_data + idx + 24, d_data + idx + 36, result2 + idx,
                   std::bind1st(std::multiplies< std::complex< double > >(), mat.d_data[2]));
    std::transform(d_data + idx,      d_data + idx + 12, result0 + idx + 12,
                   std::bind1st(std::multiplies< std::complex< double > >(), mat.d_data[3]));
    std::transform(d_data + idx + 12, d_data + idx + 24, result1 + idx + 12,
                   std::bind1st(std::multiplies< std::complex< double > >(), mat.d_data[4]));
    std::transform(d_data + idx + 24, d_data + idx + 36, result2 + idx + 12,
                   std::bind1st(std::multiplies< std::complex< double > >(), mat.d_data[5]));
    std::transform(d_data + idx,      d_data + idx + 12, result0 + idx + 24,
                   std::bind1st(std::multiplies< std::complex< double > >(), mat.d_data[6]));
    std::transform(d_data + idx + 12, d_data + idx + 24, result1 + idx + 24,
                   std::bind1st(std::multiplies< std::complex< double > >(), mat.d_data[7]));
    std::transform(d_data + idx + 24, d_data + idx + 36, result2 + idx + 24,
                   std::bind1st(std::multiplies< std::complex< double > >(), mat.d_data[8]));
  }
  // add up and write result into d_data
  std::transform(result0, result0 + 144, result1, d_data, std::plus< std::complex< double > >());
  std::transform(d_data,  d_data  + 144, result2, d_data, std::plus< std::complex< double > >());
  return *this;
}
