// #include "Tensor.ih"

inline QCD::Tensor &QCD::Tensor::rightMultiply(SU3::hcMatrix const &mat)
{
  std::complex< double > result[144];

  // this is not very efficient ...
  SU3::Matrix const mat_T = (SU3::Matrix(mat));

  std::complex< double > const zero(0, 0);

  for (size_t idx=0; idx<144; idx+=12)
  {
    result[idx +  0] = std::inner_product(mat_T.d_data,     mat_T.d_data + 3, d_data + idx,     zero);
    result[idx +  1] = std::inner_product(mat_T.d_data + 3, mat_T.d_data + 6, d_data + idx,     zero);
    result[idx +  2] = std::inner_product(mat_T.d_data + 6, mat_T.d_data + 9, d_data + idx,     zero);

    result[idx +  3] = std::inner_product(mat_T.d_data,     mat_T.d_data + 3, d_data + idx + 3, zero);
    result[idx +  4] = std::inner_product(mat_T.d_data + 3, mat_T.d_data + 6, d_data + idx + 3, zero);
    result[idx +  5] = std::inner_product(mat_T.d_data + 6, mat_T.d_data + 9, d_data + idx + 3, zero);

    result[idx +  6] = std::inner_product(mat_T.d_data,     mat_T.d_data + 3, d_data + idx + 6, zero);
    result[idx +  7] = std::inner_product(mat_T.d_data + 3, mat_T.d_data + 6, d_data + idx + 6, zero);
    result[idx +  8] = std::inner_product(mat_T.d_data + 6, mat_T.d_data + 9, d_data + idx + 6, zero);

    result[idx +  9] = std::inner_product(mat_T.d_data,     mat_T.d_data + 3, d_data + idx + 9, zero);
    result[idx + 10] = std::inner_product(mat_T.d_data + 3, mat_T.d_data + 6, d_data + idx + 9, zero);
    result[idx + 11] = std::inner_product(mat_T.d_data + 6, mat_T.d_data + 9, d_data + idx + 9, zero);
  }
  std::copy(result, result + 144, d_data);
  return *this;
}

